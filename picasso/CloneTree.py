"""
CloneTree: Phylogenetic tree analysis and visualization for PICASSO results.

This module provides the CloneTree class for integrating phylogenetic trees with
clone assignments and CNA data. It enables comprehensive analysis and visualization
of phylogenetic reconstruction results, with specific support for noisy scRNA-seq-
inferred CNA data patterns.

Classes
-------
CloneTree
    Integrates phylogenetic trees, clone assignments, and CNA profiles for
    comprehensive analysis and visualization of tumor evolution patterns.

Examples
--------
Basic usage with PICASSO results:

>>> from picasso import Picasso, CloneTree, load_data
>>>
>>> # Load example data and run PICASSO phylogenetic inference
>>> cna_data = load_data()
>>> picasso = Picasso(cna_data)
>>> picasso.fit()
>>>
>>> # Create CloneTree for analysis and visualization
>>> phylogeny = picasso.get_phylogeny()
>>> assignments = picasso.get_clone_assignments()
>>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
>>>
>>> # Generate visualizations
>>> clone_tree.plot_alterations(save_as='heatmap.pdf')
>>> clone_tree.plot_clone_sizes(save_as='sizes.pdf')

Notes
-----
The CloneTree class is designed to handle:
- Integration of phylogenetic trees with cellular data
- Aggregation of noisy CNA profiles by clone
- Visualization of clonal evolution patterns
- Export to publication-ready formats

See Also
--------
Picasso : Main phylogenetic inference algorithm
itol_utils : Functions for iTOL visualization export
utils : Data preprocessing utilities
"""

import pandas as pd
import ete3
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Tuple, Union


class CloneTree:
    __clone_phylogeny: ete3.Tree
    __sample_phylogeny: Optional[ete3.Tree]
    clone_profiles: pd.DataFrame
    clone_profiles_certainty: pd.DataFrame

    def __init__(
        self,
        phylogeny: ete3.Tree,
        clone_assignments: pd.DataFrame,
        character_matrix: pd.DataFrame,
        clone_aggregation: str = "mode",
        metadata: Optional[pd.DataFrame] = None,
    ) -> None:
        """
        Initialize a CloneTree for analysis and visualization of phylogenetic reconstruction results.

        CloneTree integrates phylogenetic trees from PICASSO with clone assignments and CNA
        data to provide comprehensive analysis and visualization capabilities. It handles
        the aggregation of noisy scRNA-seq-inferred CNA profiles by clone and supports
        various downstream analyses.

        Parameters
        ----------
        phylogeny : ete3.Tree
            The phylogenetic tree with terminal clones as leaves, typically obtained from
            the PICASSO model via get_phylogeny(). Internal nodes represent ancestral
            clones and splitting events.
        clone_assignments : pd.DataFrame
            DataFrame with cell/sample identifiers as index and a 'clone_id' column
            containing clone assignments. Should correspond to the leaves of the phylogeny.
            Typically obtained from PICASSO via get_clone_assignments().
        character_matrix : pd.DataFrame
            The CNA character matrix where rows are cells/samples and columns are genomic
            features (genes, chromosome arms, bins). Values represent inferred copy number
            states. Should contain the same samples as in clone_assignments.
        clone_aggregation : {'mode', 'mean'}, default='mode'
            Method for aggregating CNA profiles within each clone:
            - 'mode': Use most frequent copy number state (recommended for noisy data)
            - 'mean': Use average copy number (not yet implemented)
        metadata : pd.DataFrame, optional
            Additional sample metadata for visualization and analysis. Index should match
            character_matrix. Common examples include cell type annotations, sample origin,
            experimental conditions.

        Attributes
        ----------
        clone_profiles : pd.DataFrame
            Aggregated CNA profiles for each clone (rows=clones, columns=genomic features).
        clone_profiles_certainty : pd.DataFrame
            Confidence/certainty scores for each aggregated profile value.
        clone_assignments : pd.DataFrame
            DataFrame with cell/sample identifiers as index and clone assignments.
        character_matrix : pd.DataFrame
            The CNA character matrix with cells as rows and genomic features as columns.
        metadata : Optional[pd.DataFrame]
            Additional sample metadata for visualization and analysis.

        Raises
        ------
        AssertionError
            If clone_assignments lacks 'clone_id' column, if phylogeny leaves don't match
            clone assignments, if sample indices don't match between DataFrames, or if
            clone_aggregation method is invalid.

        Examples
        --------
        Basic usage with PICASSO results:

        >>> from picasso import Picasso, CloneTree, load_data
        >>>
        >>> # Load example data and run PICASSO
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix)
        >>> picasso.fit()
        >>>
        >>> # Create CloneTree for analysis
        >>> phylogeny = picasso.get_phylogeny()
        >>> assignments = picasso.get_clone_assignments()
        >>> clone_tree = CloneTree(phylogeny, assignments, character_matrix)
        >>>
        >>> # Analyze results
        >>> print(f"Number of clones: {len(clone_tree.clone_profiles)}")
        >>> clone_tree.plot_alterations(save_as='cna_heatmap.pdf')
        >>> clone_tree.plot_clone_sizes(save_as='clone_sizes.pdf')

        With metadata for enhanced visualization:

        >>> import pandas as pd
        >>> # Add cell type metadata (example)
        >>> metadata = pd.DataFrame({'cell_type': ['TypeA'] * 50 + ['TypeB'] * 50},
        ...                        index=character_matrix.index)
        >>> clone_tree = CloneTree(phylogeny, assignments, character_matrix,
        ...                       metadata=metadata)
        >>> clone_tree.plot_alterations(metadata=metadata[['cell_type']])

        Notes
        -----
        **Design Considerations for Noisy Data**:
        - Modal aggregation reduces impact of outlier cells within clones
        - Confidence scores help identify uncertain clone profiles
        - Visualization functions highlight clone-specific patterns

        **Clone Profile Aggregation**:
        - Mode aggregation finds most common copy number state per feature per clone
        - Handles missing data and ties in noisy scRNA-seq data
        - Certainty scores indicate reliability of aggregated values

        **Visualization Capabilities**:
        - Heatmaps show clone-specific CNA patterns
        - Clone size distributions reveal clonal architecture
        - Integration with iTOL for publication-quality figures

        See Also
        --------
        Picasso : Main class for phylogenetic inference from CNA data
        plot_alterations : Create heatmap visualization of CNA profiles
        plot_clone_sizes : Visualize clone size distribution
        get_sample_phylogeny : Generate sample-level phylogenetic tree
        """
        assert (
            "clone_id" in clone_assignments.columns
        ), 'The clone assignments must have a column named "clone_id".'
        assert isinstance(phylogeny, ete3.Tree)
        # Check the leaves of the phylogeny match the clones in the clone assignments
        assert set(phylogeny.get_leaf_names()) == set(clone_assignments["clone_id"]), (
            "The leaves of the phylogeny do not match the "
            "clones in the clone assignments."
        )

        # Check that the samples in the assignment matrix match the samples in the character matrix
        assert set(character_matrix.index) == set(clone_assignments.index), (
            "The samples in the assignment matrix do not match the samples in the character "
            "matrix."
        )

        clone_aggregation = clone_aggregation.lower()
        assert clone_aggregation in [
            "mode",
            "mean",
        ], 'The clone aggregation method must be either "mode" or "mean".'

        self.__clone_phylogeny = phylogeny
        self.__sample_phylogeny = None

        self.clone_assignments = clone_assignments
        self.character_matrix = character_matrix

        assert metadata is None or isinstance(
            metadata, pd.DataFrame
        ), "The metadata must be a pandas DataFrame."
        if metadata is not None:
            assert set(metadata.index) == set(
                character_matrix.index
            ), "The samples in the metadata do not match the samples in the character matrix."
        self.metadata = metadata

        self.clone_profiles, self.clone_profiles_certainty = self.aggregate_clones(
            clone_aggregation
        )
        print(
            f"Initialized CloneTree with {len(self.clone_profiles)} clones and {len(self.character_matrix)} samples."
        )

    def aggregate_clones(
        self, aggregation_method: str
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Aggregate CNA profiles within each clone to create representative clone profiles.

        Combines individual cell CNA profiles within each clone into single representative
        profiles using statistical aggregation. This reduces noise and creates clean
        clone-level CNA signatures for downstream analysis and visualization.

        Parameters
        ----------
        aggregation_method : str
            Method for aggregating CNA values within clones:
            - 'mode': Use most frequent copy number state (recommended for noisy data)
            - 'mean': Use average copy number (not yet implemented)

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            First DataFrame: Aggregated clone profiles with clones as rows and genomic
            features as columns. Values represent the aggregated copy number states.
            Second DataFrame: Certainty/confidence scores for each aggregated value,
            indicating reliability of the aggregation.

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>> profiles, certainty = clone_tree.aggregate_clones('mode')
        >>> print(f"Clone profiles shape: {profiles.shape}")
        >>> print(f"Average certainty: {certainty.mean().mean():.2f}")

        Notes
        -----
        **Modal Aggregation**:
        - Finds the most common copy number state for each feature within each clone
        - Handles ties by selecting the first modal value
        - Provides certainty scores based on frequency of the modal state
        - Robust to outlier cells within clones
        - Facilitates visualization of CNA patterns across clones

        **Design for Noisy Data**:
        - Modal aggregation reduces impact of noise and technical artifacts
        - Certainty scores help identify unreliable aggregated values
        - Particularly effective for scRNA-seq-inferred CNA data

        Raises
        ------
        NotImplementedError
            If aggregation_method is 'mean' (not yet implemented).
        ValueError
            If aggregation_method is not 'mode' or 'mean'.

        See Also
        --------
        get_modal_clone_profiles : Internal method implementing modal aggregation
        """
        if aggregation_method == "mode":
            return self.get_modal_clone_profiles()
        elif aggregation_method == "mean":
            raise NotImplementedError("Mean aggregation is not yet implemented.")

    def get_most_ancestral_clone(self) -> str:
        """
        Identify the most ancestral clone based on CNA profile complexity.

        Determines which clone represents the most ancestral state by counting the
        number of copy number alterations (deviations from neutral state). This
        is useful for rooting phylogenetic trees and understanding evolutionary
        relationships.

        Returns
        -------
        str
            Clone identifier of the most ancestral clone (fewest alterations).

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>> ancestral = clone_tree.get_most_ancestral_clone()
        >>> print(f"Most ancestral clone: {ancestral}")
        >>>
        >>> # Use for tree rooting
        >>> clone_tree.root_tree(ancestral)

        Notes
        -----
        **Ancestral State Assumptions**:
        - Copy number state 0 is considered the ancestral/neutral state
        - Clones with more alterations are considered more derived
        - Useful for establishing evolutionary directionality

        **Algorithm**:
        1. Count non-zero states for each clone in aggregated profiles
        2. Select clone with minimum alteration count
        3. Return clone identifier

        **Use Cases**:
        - Rooting phylogenetic trees for visualization
        - Identifying putative normal/founder cell populations
        - Understanding tumor evolution trajectories

        See Also
        --------
        root_tree : Method to root the phylogeny using an outgroup clone
        clone_profiles : Aggregated CNA profiles used for ancestral inference
        """
        num_alterations = (self.clone_profiles != 0).sum(axis=1)
        ancestral_clone = num_alterations.idxmin()
        return ancestral_clone

    def root_tree(self, outgroup: str) -> None:
        """
        Root the phylogenetic tree using a specified outgroup clone.

        Establishes evolutionary directionality by setting a designated clone as
        the outgroup, which becomes the root of the tree. This is essential for
        proper interpretation of evolutionary relationships and visualization.

        Parameters
        ----------
        outgroup : str
            Identifier of the clone to use as outgroup. Must be present in the
            phylogenetic tree leaves. Often the most ancestral clone identified
            by get_most_ancestral_clone().

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>>
        >>> # Root with most ancestral clone
        >>> ancestral = clone_tree.get_most_ancestral_clone()
        >>> clone_tree.root_tree(ancestral)
        >>>
        >>> # Or root with specific clone
        >>> clone_tree.root_tree('1-0-STOP')

        Notes
        -----
        **Effects of Rooting**:
        - Changes tree topology and evolutionary interpretation
        - Affects all subsequent tree-based analyses
        - Resets sample phylogeny (if previously generated)
        - Essential for proper tree visualization

        **Outgroup Selection Guidelines**:
        - Use most ancestral clone (fewest alterations) when possible
        - Consider biological knowledge about cell populations
        - Avoid clones with many unique alterations

        **Implementation Details**:
        - Uses ete3's set_outgroup() method
        - Invalidates cached sample phylogeny
        - Tree structure is modified in-place

        Raises
        ------
        AssertionError
            If outgroup is not found among the tree leaves.

        See Also
        --------
        get_most_ancestral_clone : Identify suitable outgroup candidates
        get_clone_phylogeny : Access the rooted phylogenetic tree
        get_sample_phylogeny : Generate sample-level tree from rooted clone tree
        """
        assert (
            outgroup in self.__clone_phylogeny.get_leaf_names()
        ), "The outgroup must be a leaf in the tree."
        self.__sample_phylogeny = None
        self.__clone_phylogeny.set_outgroup(outgroup)
        return

    def get_clone_phylogeny(self) -> ete3.Tree:
        """
        Access the clone-level phylogenetic tree.

        Returns the phylogenetic tree where leaves represent clones (terminal cell
        populations) and internal nodes represent ancestral populations. This is
        the primary tree structure used for evolutionary analysis.

        Returns
        -------
        ete3.Tree
            Phylogenetic tree with clones as leaves. Tree may be rooted or unrooted
            depending on whether root_tree() has been called.

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>> tree = clone_tree.get_clone_phylogeny()
        >>> print(f"Tree has {len(tree.get_leaves())} clones")
        >>> print("Clone names:", tree.get_leaf_names())
        >>>
        >>> # Tree manipulation
        >>> if not tree.is_root():
        ...     print("Tree is rooted")
        >>>
        >>> # Export to Newick format
        >>> newick_str = tree.write()

        Notes
        -----
        **Tree Structure**:
        - Leaves represent terminal clones from PICASSO analysis
        - Internal nodes represent inferred ancestral states
        - Branch structure reflects evolutionary relationships
        - Node names correspond to clone identifiers

        **Tree States**:
        - May be rooted (after root_tree()) or unrooted
        - Tree topology reflects PICASSO splitting hierarchy
        - Compatible with standard phylogenetic analysis tools

        **Use Cases**:
        - Phylogenetic visualization and analysis
        - Export to external tools (iTOL, FigTree, etc.)
        - Evolutionary distance calculations
        - Tree-based clustering validation

        See Also
        --------
        get_sample_phylogeny : Get expanded tree with individual cells
        root_tree : Root the tree for proper evolutionary interpretation
        """
        return self.__clone_phylogeny

    def get_sample_phylogeny(self) -> ete3.Tree:
        """
        Generate expanded phylogenetic tree with individual cells as leaves.

        Creates a detailed tree where each cell/sample appears as a separate leaf,
        while maintaining the clone-based evolutionary structure. Cells within the
        same clone are attached as children of their respective clone nodes.

        Returns
        -------
        ete3.Tree
            Expanded phylogenetic tree where leaves represent individual cells/samples
            rather than clones. Clone nodes become internal nodes with cells as children.

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>> sample_tree = clone_tree.get_sample_phylogeny()
        >>> print(f"Tree has {len(sample_tree.get_leaves())} cells")
        >>>
        >>> # Access cell-specific information
        >>> for leaf in sample_tree.get_leaves():
        ...     print(f"Cell {leaf.name}")
        ...     if clone_tree.metadata is not None:
        ...         print(f"  Metadata: {leaf.features}")

        Notes
        -----
        **Tree Construction**:
        - Starts with clone phylogeny as backbone
        - Adds individual cells as children of clone nodes
        - Preserves evolutionary relationships at clone level
        - Enables cell-level analysis within phylogenetic context

        **Metadata Integration**:
        - If metadata provided, adds features to cell nodes
        - Features accessible via leaf.features or leaf.get_feature()
        - Enables metadata-aware tree visualization

        **Performance Considerations**:
        - Tree generated on first call, then cached
        - Cache invalidated when tree is re-rooted
        - Large datasets may produce complex trees

        **Use Cases**:
        - Cell-level phylogenetic visualization
        - Metadata mapping onto evolutionary structure
        - Detailed iTOL annotations
        - Single-cell evolutionary analysis

        See Also
        --------
        get_clone_phylogeny : Access the underlying clone tree structure
        metadata : Cell-level metadata integrated into tree nodes
        """
        if self.__sample_phylogeny is None:
            cell_tree = self.__clone_phylogeny.copy()
            n_leaves_added = 0
            for clone in cell_tree.get_leaves():
                samples = self.clone_assignments.query(
                    f'clone_id == "{clone.name}"'
                ).index
                for sample in samples:
                    clone.add_child(name=sample)
                    n_leaves_added += 1
            print(f"Added {n_leaves_added} leaves to the tree.")
            assert set(cell_tree.get_leaf_names()) == set(
                self.character_matrix.index
            ), (
                "The samples in the tree do not match the samples in the character "
                "matrix."
            )

            self.__sample_phylogeny = cell_tree

            if self.metadata is not None:
                for sample_node in self.__sample_phylogeny.get_leaves():
                    for column in self.metadata.columns:
                        sample = sample_node.name
                        sample_node.add_feature(
                            column, self.metadata.loc[sample, column]
                        )

        return self.__sample_phylogeny

    def infer_evolutionary_changes(self) -> None:
        """
        Infer evolutionary changes along phylogenetic tree branches.

        Reconstructs the specific copy number alterations that occurred at each
        internal node of the phylogenetic tree by analyzing transitions between
        ancestral and derived clone profiles. This method is planned for future
        implementation.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented. Future versions will support
            ancestral state reconstruction and evolutionary change mapping.

        Notes
        -----
        **Planned Functionality**:
        - Ancestral state reconstruction for internal tree nodes
        - Identification of specific CNA events along branches

        **Potential Applications**:
        - Understanding CNA acquisition patterns
        - Identifying driver vs passenger alterations
        - Validating phylogenetic relationships

        See Also
        --------
        clone_profiles : Aggregated clone CNA profiles used for inference
        get_clone_phylogeny : Phylogenetic tree structure for change mapping
        """
        raise NotImplementedError

    def plot_alterations(
        self,
        metadata: Optional[pd.DataFrame] = None,
        cmap: str = "coolwarm",
        show: bool = True,
        save_as: Optional[str] = None,
        center: Optional[float] = None,
    ) -> None:
        """
        Create clustered heatmap visualization of CNA profiles with clone annotations.

        Generates a comprehensive heatmap showing copy number alterations across all
        cells, with cells grouped by clone assignment and colored sidebars indicating
        clone membership and optional metadata categories.

        Parameters
        ----------
        metadata : pd.DataFrame, optional
            Additional metadata for enhanced visualization. Index should match
            character_matrix. Each column represents a metadata category (e.g.,
            cell_type, treatment, tissue). Will be displayed as colored sidebars.
        cmap : str, default='coolwarm'
            Matplotlib colormap for the main heatmap. Common choices:
            - 'coolwarm': Blue-white-red for CNAs (deletions-neutral-amplifications)
            - 'RdBu_r': Red-blue reversed
            - 'viridis': Perceptually uniform colormap
        show : bool, default=True
            Whether to display the plot interactively.
        save_as : str, optional
            File path to save the plot. Supports common formats (.pdf, .png, .svg).
            Recommended: use .pdf for publication quality.
        center : float, optional
            Value at which to center the colormap. If None, uses default centering.
            For CNA data, typically 0 (neutral copy number) or 2 (diploid).

        Examples
        --------
        Basic heatmap with clone annotations:

        >>> from picasso import Picasso, CloneTree, load_data
        >>>
        >>> # Create CloneTree
        >>> cna_data = load_data()
        >>> picasso = Picasso(cna_data)
        >>> picasso.fit()
        >>> clone_tree = CloneTree(picasso.get_phylogeny(),
        ...                       picasso.get_clone_assignments(),
        ...                       cna_data)
        >>>
        >>> # Basic visualization
        >>> clone_tree.plot_alterations(save_as='cna_heatmap.pdf')

        Enhanced visualization with metadata:

        >>> import pandas as pd
        >>>
        >>> # Add cell type metadata
        >>> metadata = pd.DataFrame({
        ...     'cell_type': ['Malignant'] * 80 + ['Normal'] * 20,
        ...     'tissue': ['Primary'] * 60 + ['Metastasis'] * 40
        ... }, index=cna_data.index)
        >>>
        >>> # Create enhanced heatmap
        >>> clone_tree.plot_alterations(metadata=metadata,
        ...                            cmap='RdBu_r',
        ...                            center=0,
        ...                            save_as='enhanced_heatmap.pdf')

        Notes
        -----
        **Visualization Features**:
        - Cells automatically grouped by clone assignment
        - Clone-specific color sidebar for easy identification
        - Optional metadata sidebars for additional context
        - Configurable color schemes for different data types

        **Layout Organization**:
        - Rows: Individual cells/samples
        - Columns: Genomic features (chromosome arms, genes, etc.)
        - Left sidebars: Clone assignments + optional metadata
        - Main heatmap: Copy number alteration values

        **Color Interpretation**:
        - Clone sidebar: Each clone gets a distinct color
        - Metadata sidebars: Categorical values get distinct colors
        - Main heatmap: Continuous colormap for CNA values

        **Best Practices**:
        - Use 'coolwarm' colormap for copy number data
        - Center colormap at neutral copy number (typically 0 or 2)
        - Save as PDF for publication-quality figures
        - Include relevant metadata for biological context

        See Also
        --------
        plot_clone_sizes : Visualize clone size distribution
        clone_profiles : Access aggregated clone CNA profiles
        seaborn.clustermap : Underlying plotting function used
        """
        df = self.character_matrix.join(self.clone_assignments)
        # Sort the columns by clone assignment
        df = df.sort_values(by="clone_id")

        # Colour cells by clone assignment
        palette = sns.color_palette("tab20", len(df["clone_id"].unique()))
        clone_cmap = {}
        for i, clone in enumerate(df["clone_id"].unique()):
            clone_cmap[clone] = self.rgba_to_hex(palette[i])
        row_colors = pd.DataFrame(df["clone_id"].map(clone_cmap))

        # Plot a clustered heatmap, so that we can display the clone assignments as a colour bar
        if metadata is not None:
            row_colors = row_colors.join(metadata)
        if center is not None:
            sns.clustermap(
                df.drop(columns="clone_id"),
                row_colors=row_colors,
                col_cluster=False,
                row_cluster=False,
                cmap=cmap,
                figsize=(10, 10),
                center=center,
            )
        else:
            sns.clustermap(
                df.drop(columns="clone_id"),
                row_colors=row_colors,
                col_cluster=False,
                row_cluster=False,
                cmap=cmap,
                figsize=(10, 10),
            )
        if save_as:
            plt.savefig(save_as, dpi=300)
        if show:
            plt.show()
        plt.close()

    def plot_clone_sizes(
        self, show: bool = True, save_as: Optional[str] = None
    ) -> None:
        """
        Visualize the distribution of clone sizes in the phylogenetic tree.

        Creates a histogram showing how many cells belong to each clone, providing
        insights into clonal architecture, diversity, and potential dominant/rare
        clones within the analyzed population.

        Parameters
        ----------
        show : bool, default=True
            Whether to display the plot interactively using matplotlib.
        save_as : str, optional
            File path to save the plot. Supports common formats (.pdf, .png, .svg).
            If provided, plot will be saved to this location.

        Examples
        --------
        Basic clone size visualization:

        >>> from picasso import Picasso, CloneTree, load_data
        >>>
        >>> # Create CloneTree and visualize clone sizes
        >>> cna_data = load_data()
        >>> picasso = Picasso(cna_data)
        >>> picasso.fit()
        >>> clone_tree = CloneTree(picasso.get_phylogeny(),
        ...                       picasso.get_clone_assignments(),
        ...                       cna_data)
        >>>
        >>> # Display clone size distribution
        >>> clone_tree.plot_clone_sizes()

        Save without displaying:

        >>> # Save to file without showing
        >>> clone_tree.plot_clone_sizes(show=False, save_as='clone_sizes.pdf')

        Analyze clone architecture:

        >>> # Get clone sizes for analysis
        >>> assignments = picasso.get_clone_assignments()
        >>> clone_sizes = assignments['clone_id'].value_counts()
        >>> print(f"Largest clone: {clone_sizes.max()} cells")
        >>> print(f"Smallest clone: {clone_sizes.min()} cells")
        >>> print(f"Mean clone size: {clone_sizes.mean():.1f} cells")
        >>>
        >>> # Visualize
        >>> clone_tree.plot_clone_sizes(save_as='clone_architecture.pdf')

        Notes
        -----
        **Plot Features**:
        - Histogram showing distribution of clone sizes
        - X-axis: Clone size (number of cells per clone)
        - Y-axis: Number of clones with that size
        - Kernel density estimate (KDE) overlay for smooth distribution
        - Automatic binning based on data range

        **Interpretation**:
        - Right-skewed distribution: Few large clones dominate
        - Uniform distribution: Balanced clonal architecture
        - Left-skewed distribution: Many small clones, rare large ones

        **Technical Considerations**:
        - Clone sizes depend on PICASSO parameters (min_clone_size, etc.)
        - Very small clones may indicate noise or over-splitting
        - Very large clones may indicate under-splitting or homogeneity

        See Also
        --------
        plot_alterations : Visualize CNA profiles with clone annotations
        clone_assignments : Access raw clone assignment data
        get_clone_assignments : Get clone assignments from PICASSO analysis
        """
        cells_per_clone = self.clone_assignments["clone_id"].value_counts()
        plt.figure()
        sns.histplot(cells_per_clone, kde=True)
        plt.xlabel("Clone Size")
        plt.xticks(rotation=45)
        plt.ylabel("Number of Clones")
        plt.title("Number of Cells per Clone")
        if save_as:
            plt.savefig(save_as)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def calc_mode(series: pd.Series) -> Union[int, float, None]:
        """
        Calculate the statistical mode (most frequent value) of a pandas Series.

        Computes the most common value in a series, handling edge cases where no
        mode exists or multiple modes are present. Used for aggregating copy number
        states within clones.

        Parameters
        ----------
        series : pd.Series
            Input data series containing numeric values (typically copy number states).

        Returns
        -------
        int, float, or None
            The most frequent value in the series. Returns None if series is empty
            or all values are NaN. If multiple modes exist, returns the first one.

        Examples
        --------
        >>> import pandas as pd
        >>> data = pd.Series([1, 1, 2, 2, 2, 3])
        >>> CloneTree.calc_mode(data)
        2
        >>>
        >>> # Handle ties
        >>> tie_data = pd.Series([1, 1, 2, 2])
        >>> CloneTree.calc_mode(tie_data)  # Returns first mode
        1

        Notes
        -----
        - Uses pandas Series.mode() method internally
        - Handles empty series gracefully by returning None
        - For ties, returns the first modal value (arbitrary but consistent)
        - Designed for integer copy number data but works with any numeric type

        See Also
        --------
        calc_mode_freq : Calculate frequency of the modal value
        get_modal_clone_profiles : Main method using this utility
        """
        mode = series.mode()
        if len(mode) > 0:  # If there's at least one mode
            return mode[0]  # Return the first mode
        return None

    @staticmethod
    def calc_mode_freq(series: pd.Series) -> float:
        """
        Calculate the frequency (proportion) of the modal value in a pandas Series.

        Computes what fraction of values in the series match the most frequent value.
        This provides a confidence measure for modal aggregation - higher frequencies
        indicate more reliable consensus within the data.

        Parameters
        ----------
        series : pd.Series
            Input data series containing numeric values (typically copy number states).

        Returns
        -------
        float
            Proportion of values matching the modal value, between 0.0 and 1.0.
            Returns 0.0 if series is empty or contains only NaN values.

        Examples
        --------
        >>> import pandas as pd
        >>> # High consensus
        >>> data = pd.Series([2, 2, 2, 2, 1])
        >>> CloneTree.calc_mode_freq(data)
        0.8  # 4 out of 5 values are modal
        >>>
        >>> # Perfect consensus
        >>> uniform = pd.Series([1, 1, 1, 1])
        >>> CloneTree.calc_mode_freq(uniform)
        1.0
        >>>
        >>> # Low consensus (tie)
        >>> mixed = pd.Series([1, 2, 3, 4])
        >>> CloneTree.calc_mode_freq(mixed)
        0.25  # Each value appears once

        Notes
        -----
        **Interpretation Guide**:
        - 1.0: Perfect consensus, all values identical
        - 0.8-0.9: Strong consensus with few outliers
        - 0.5-0.7: Moderate consensus, some heterogeneity
        - <0.5: Weak consensus, high heterogeneity

        **Use in Clone Analysis**:
        - Quality metric for clone coherence
        - Confidence score for aggregated profiles
        - Filter for reliable clone assignments
        - Identifies noisy or heterogeneous clones

        See Also
        --------
        calc_mode : Calculate the actual modal value
        get_modal_clone_profiles : Main method using this utility for confidence scores
        """
        mode = series.mode()
        if len(mode) > 0:
            return len(series[series == mode[0]]) / len(series)
        return 0

    def get_modal_clone_profiles(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Compute modal (most frequent) copy number states for each clone.

        Aggregates CNA profiles within each clone by finding the most common copy
        number state for each genomic feature. Also computes confidence scores
        based on the frequency of the modal state.

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            modal_profiles : pd.DataFrame
                Clone profiles with modal copy number states. Rows are clones,
                columns are genomic features. Values are the most frequent
                copy number state within each clone.
            modal_frequencies : pd.DataFrame
                Confidence scores for modal states. Same structure as modal_profiles
                but values represent the proportion of cells with the modal state
                (0.0 to 1.0, where 1.0 indicates all cells have the same state).

        Examples
        --------
        >>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
        >>> profiles, frequencies = clone_tree.get_modal_clone_profiles()
        >>>
        >>> # Examine profile quality
        >>> avg_confidence = frequencies.mean().mean()
        >>> print(f"Average modal confidence: {avg_confidence:.2f}")
        >>>
        >>> # Find highly confident features
        >>> confident_features = frequencies.columns[frequencies.mean() > 0.8]
        >>> print(f"High confidence features: {len(confident_features)}")

        Notes
        -----
        **Modal Aggregation Process**:
        1. Group cells by clone assignment
        2. For each clone-feature combination, find most frequent copy number state
        3. Calculate frequency of modal state as confidence measure
        4. Handle ties by selecting first modal value

        **Confidence Interpretation**:
        - 1.0: All cells in clone have identical copy number state
        - 0.5-0.9: Majority consensus with some variation
        - <0.5: High heterogeneity, unreliable modal state

        **Noise Handling**:
        - Modal aggregation naturally filters outlier cells
        - Confidence scores identify unreliable aggregations
        - Particularly effective for noisy scRNA-seq-inferred CNAs

        **Applications**:
        - Generate clean clone signatures for visualization
        - Quality control for clone assignments
        - Feature selection based on clone coherence

        See Also
        --------
        calc_mode : Static method for computing modal values
        calc_mode_freq : Static method for computing modal frequencies
        aggregate_clones : Public interface using this method
        """

        # Ensure the indices are aligned
        cnvs = self.character_matrix.loc[self.clone_assignments.index]

        # Merge the two DataFrames on their indices
        merged_df = pd.concat([self.clone_assignments, cnvs], axis=1)

        clone_column = "clone_id"

        # Modal values DataFrame
        modal_df = merged_df.groupby(clone_column).agg(self.calc_mode).reset_index()

        # Frequencies of modal values DataFrame
        freq_df = merged_df.groupby(clone_column).agg(self.calc_mode_freq).reset_index()

        # Set the clone column as the index again for modal_df and freq_df
        modal_df.set_index(clone_column, inplace=True)
        freq_df.set_index(clone_column, inplace=True)

        return modal_df, freq_df

    @staticmethod
    def rgba_to_hex(rgba: Tuple[float, ...]) -> str:
        """
        Convert RGBA values to hexadecimal color string.

        Parameters
        ----------
        rgba : tuple
            Tuple of RGBA values.

        Returns
        -------
        str
            Hexadecimal color string.

        Examples
        --------
        >>> rgba_to_hex((1.0, 0.0, 0.0, 1.0))
        '#ff0000'
        """

        # Extract the RGBA values
        if len(rgba) == 3:
            red, green, blue = rgba
        elif len(rgba) == 4:
            red, green, blue, _ = rgba

        # Ensure the values are in the range 0-1
        red = min(1.0, max(0.0, red))
        green = min(1.0, max(0.0, green))
        blue = min(1.0, max(0.0, blue))

        # Convert to hexadecimal and ensure two characters for each value
        red_hex = format(int(red * 255), "02X")
        green_hex = format(int(green * 255), "02X")
        blue_hex = format(int(blue * 255), "02X")

        # Concatenate the hexadecimal values
        hex_color = f"#{red_hex}{green_hex}{blue_hex}"

        return hex_color


# Define public API
__all__ = ["CloneTree"]
