import pandas as pd
import ete3
import seaborn as sns
import matplotlib.pyplot as plt


class CloneTree:
    def __init__(self, phylogeny, clone_assignments, character_matrix, clone_aggregation='mode', metadata=None):
        """
        Initialize a CloneTree object. This object represents a phylogenetic tree with clones as leaves and allows for
        the aggregation of character matrices by clone.

        :param phylogeny: (ete3.Tree) The phylogenetic tree with clones as leaves.
        :param clone_assigments: (pandas.DataFrame) A DataFrame with the clone assignments for each sample. The index
        should be the sample names and there should be a column named 'clone_id' with the clone assignments.
        :param character_matrix: (pandas.DataFrame) A DataFrame with the character matrix. Each row should be a sample
        and each column should be a feature representing some integer genomic alteration.
        :param clone_aggregation: (str) The method to use for aggregating the character matrix by clone. Either 'mode'
        or 'mean'.
        :param metadata: (pandas.DataFrame) Optional metadata for the samples.
        """
        assert 'clone_id' in clone_assignments.columns, 'The clone assignments must have a column named "clone_id".'
        assert isinstance(phylogeny, ete3.Tree)
        # Check the leaves of the phylogeny match the clones in the clone assignments
        assert set(phylogeny.get_leaf_names()) == set(clone_assignments[
                                                          'clone_id']), ('The leaves of the phylogeny do not match the '
                                                                         'clones in the clone assignments.')

        # Check that the samples in the assignment matrix match the samples in the character matrix
        assert set(character_matrix.index) == set(
            clone_assignments.index), ('The samples in the assignment matrix do not match the samples in the character '
                                       'matrix.')

        clone_aggregation = clone_aggregation.lower()
        assert clone_aggregation in ['mode', 'mean'], 'The clone aggregation method must be either "mode" or "mean".'

        self.__clone_phylogeny = phylogeny
        self.__sample_phylogeny = None

        self.clone_assignments = clone_assignments
        self.character_matrix = character_matrix

        assert metadata is None or isinstance(metadata, pd.DataFrame), 'The metadata must be a pandas DataFrame.'
        if metadata is not None:
            assert set(metadata.index) == set(
                character_matrix.index), 'The samples in the metadata do not match the samples in the character matrix.'
        self.metadata = metadata

        self.clone_profiles, self.clone_profiles_certainty = self.aggregate_clones(clone_aggregation)
        print(f'Initialized CloneTree with {len(self.clone_profiles)} clones and {len(self.character_matrix)} samples.')

    def aggregate_clones(self, aggregation_method):
        """
        Aggregate the character matrix by clone to get a clone profile matrix
        :param aggregation_method: (str) The method to use for aggregation. Either 'mode' or 'mean'.
        :return:
        """
        if aggregation_method == 'mode':
            return self.get_modal_clone_profiles()
        elif aggregation_method == 'mean':
            raise NotImplementedError('Mean aggregation is not yet implemented.')

    def get_most_ancestral_clone(self):
        """
        Get the most ancestral clone in the clone tree. This is the node with the fewest number of alterations in the
        character matrix. We assume the 0 state is the ancestral state.

        :return: (str) The name of the most ancestral clone.
        """
        num_alterations = (self.clone_profiles != 0).sum(axis=1)
        ancestral_clone = num_alterations.idxmin()
        return ancestral_clone

    def root_tree(self, outgroup):
        """
        Root the tree with the given outgroup.

        :param outgroup: (str) The name of the outgroup sample.
        :return: None
        """
        assert outgroup in self.__clone_phylogeny.get_leaf_names(), 'The outgroup must be a leaf in the tree.'
        self.__sample_phylogeny = None
        self.__clone_phylogeny.set_outgroup(outgroup)
        return

    def get_clone_phylogeny(self):
        """
        Get the tree with clones as leaves.
        :return:
        """
        return self.__clone_phylogeny

    def get_sample_phylogeny(self):
        """
        Get the tree with samples as leaves.
        :return:
        """
        if self.__sample_phylogeny is None:
            cell_tree = self.__clone_phylogeny.copy()
            n_leaves_added = 0
            for clone in cell_tree.get_leaves():
                samples = self.clone_assignments.query(f'clone_id == "{clone.name}"').index
                for sample in samples:
                    clone.add_child(name=sample)
                    n_leaves_added += 1
            print(f'Added {n_leaves_added} leaves to the tree.')
            assert set(cell_tree.get_leaf_names()) == set(
                    self.character_matrix.index), ('The samples in the tree do not match the samples in the character '
                                                   'matrix.')

            self.__sample_phylogeny = cell_tree

            if self.metadata is not None:
                for sample_node in self.__sample_phylogeny.get_leaves():
                    for column in self.metadata.columns:
                        sample = sample_node.name
                        sample_node.add_feature(column, self.metadata.loc[sample, column])

        return self.__sample_phylogeny

    def infer_evolutionary_changes(self):
        """
        Infer the evolutionary changes that occurred on the tree.

        :return:
        """
        raise NotImplementedError

    def plot_alterations(self, metadata=None, cmap='coolwarm', show=True, save_as=None, center=None):
        """
        Plot the alterations in a heatmap, coloured by clone assignment and other potential metadata. :param
        color_metadata: (pandas.DataFrame) The metadata to colour the heatmap by. This should have the same index as
        the character matrix. Each column is a category and the values are the colours for each sample. :param show:
        (bool) Whether to show the plot. :param save_as: (str) The path to save the plot to. :return:
        """
        df = self.character_matrix.join(self.clone_assignments)
        # Sort the columns by clone assignment
        df = df.sort_values(by='clone_id')

        # Colour cells by clone assignment
        palette = sns.color_palette('tab20', len(df['clone_id'].unique()))
        clone_cmap = {}
        for i, clone in enumerate(df['clone_id'].unique()):
            clone_cmap[clone] = self.rgba_to_hex(palette[i])
        row_colors = pd.DataFrame(df['clone_id'].map(clone_cmap))

        # Plot a clustered heatmap, so that we can display the clone assignments as a colour bar
        if metadata is not None:
            row_colors = row_colors.join(metadata)
        if center is not None:
            sns.clustermap(df.drop(columns='clone_id'), row_colors=row_colors, col_cluster=False, row_cluster=False,
                           cmap=cmap, figsize=(10, 10), center=center)
        else:
            sns.clustermap(df.drop(columns='clone_id'), row_colors=row_colors, col_cluster=False, row_cluster=False,
                           cmap=cmap, figsize=(10, 10))
        if save_as:
            plt.savefig(save_as, dpi=300)
        if show:
            plt.show()
        plt.close()

    def plot_clone_sizes(self, show=True, save_as=None):
        """
        Plot the sizes of the clones in the tree.
        :param show: (bool) Whether to show the plot.
        :param save_as: (str) The path to save the plot to.
        :return:
        """
        cells_per_clone = self.clone_assignments['clone_id'].value_counts()
        plt.figure()
        sns.histplot(cells_per_clone, kde=True)
        plt.xlabel('Clone Size')
        plt.xticks(rotation=45)
        plt.ylabel('Number of Clones')
        plt.title('Number of Cells per Clone')
        if save_as:
            plt.savefig(save_as)
        if show:
            plt.show()
        plt.close()

    @staticmethod
    def calc_mode(series):
        mode = series.mode()
        if len(mode) > 0:  # If there's at least one mode
            return mode[0]  # Return the first mode
        return None

    @staticmethod
    def calc_mode_freq(series):
        mode = series.mode()
        if len(mode) > 0:
            return len(series[series == mode[0]]) / len(series)
        return 0

    def get_modal_clone_profiles(self):
        """
        Get the modal clone profiles for each cell. :param character_matrix: pandas.DataFrame, where each row is a
        cell and each column is a feature representing some integer genomic alteration :param clones:
        pandas.DataFrame, where each row is a cell and each column is its clone assignment :return: pandas.DataFrame,
        where each row is a clone and each column is a genomic region containing the modal copy number
        pandas.DataFrame, where each row is a clone and each column is a genomic region containing the frequency of
        the modal copy number
        """

        # Ensure the indices are aligned
        cnvs = self.character_matrix.loc[self.clone_assignments.index]

        # Merge the two DataFrames on their indices
        merged_df = pd.concat([self.clone_assignments, cnvs], axis=1)

        clone_column = 'clone_id'

        # Modal values DataFrame
        modal_df = merged_df.groupby(clone_column).agg(self.calc_mode).reset_index()

        # Frequencies of modal values DataFrame
        freq_df = merged_df.groupby(clone_column).agg(self.calc_mode_freq).reset_index()

        # Set the clone column as the index again for modal_df and freq_df
        modal_df.set_index(clone_column, inplace=True)
        freq_df.set_index(clone_column, inplace=True)

        return modal_df, freq_df

    @staticmethod
    def rgba_to_hex(rgba):
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
        red_hex = format(int(red * 255), '02X')
        green_hex = format(int(green * 255), '02X')
        blue_hex = format(int(blue * 255), '02X')

        # Concatenate the hexadecimal values
        hex_color = f"#{red_hex}{green_hex}{blue_hex}"

        return hex_color
