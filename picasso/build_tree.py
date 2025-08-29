"""
PICASSO: Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization.

This module implements the core PICASSO algorithm for reconstructing tumor phylogenies
from noisy, inferred copy number alteration (CNA) data derived from single-cell RNA
sequencing. The algorithm uses iterative binary splitting with categorical mixture
models to handle uncertainty and noise typical in scRNA-seq-inferred CNAs.

Classes
-------
Picasso
    Main class implementing the phylogenetic inference algorithm with noise handling
    capabilities designed specifically for scRNA-seq-inferred CNA data.

Examples
--------
Basic phylogenetic reconstruction:

>>> from picasso import Picasso, load_data
>>>
>>> # Load example CNA data
>>> cna_data = load_data()
>>>
>>> # Initialize with parameters suitable for noisy data
>>> picasso = Picasso(cna_data,
...                  min_clone_size=10,  # Larger for noisy data
...                  assignment_confidence_threshold=0.8)
>>>
>>> # Reconstruct phylogeny
>>> picasso.fit()
>>> phylogeny = picasso.get_phylogeny()
>>> assignments = picasso.get_clone_assignments()

Notes
-----
The PICASSO algorithm is specifically designed to handle the challenges of:
- Noise and artifacts in scRNA-seq-inferred CNAs
- Uncertainty in copy number state assignments
- Variable clone sizes and imbalanced data
- Over-fitting to noise patterns

See Also
--------
CloneTree : Visualization and analysis of phylogenetic results
utils : Utility functions for data preprocessing and loading
itol_utils : Functions for creating iTOL-compatible visualizations
"""

from pomegranate.gmm import GeneralMixtureModel
from pomegranate.distributions import Categorical

import copy
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
from tqdm import tqdm
import time

import ete3

# Set up logging
import logging

log = logging.getLogger()
log.setLevel(logging.INFO)
if not log.hasHandlers():
    log.addHandler(logging.StreamHandler())
log.propagate = False


class Picasso:
    terminal_clones: Dict[str, pd.Index]
    clones: Dict[str, pd.Index]
    depth: int

    def __init__(
        self,
        character_matrix: pd.DataFrame,
        min_depth: Optional[int] = None,
        max_depth: Optional[int] = None,
        min_clone_size: int = 5,
        terminate_by: str = "probability",
        assignment_confidence_threshold: float = 0.75,
        assignment_confidence_proportion: float = 0.8,
        bic_penalty_strength: float = 1.0,
    ) -> None:
        """
        Initialize the PICASSO model for phylogenetic inference from noisy CNA data.

        PICASSO (Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data
        Optimization) reconstructs tumor phylogenies from inferred copy number alterations (CNAs)
        derived from single-cell RNA sequencing data. Unlike direct scDNA-seq data, scRNA-seq-inferred
        CNAs are noisy and require specialized handling for more accurate phylogenetic reconstruction.

        Parameters
        ----------
        character_matrix : pd.DataFrame
            An integer matrix where rows are single cells/samples and columns are genomic features
            (e.g., chromosome arms, genes, or genomic bins). Values represent inferred copy number
            states (e.g., 0=deletion, 1=neutral, 2=amplification). For noisy scRNA-seq-inferred
            data, values may include noise artifacts that PICASSO handles through probabilistic
            modeling.
        min_depth : int, optional
            The minimum depth (number of splitting iterations) of the phylogeny. Forces algorithm
            to continue splitting even if termination criteria are met, useful for exploring
            deeper clonal structure in noisy data. Default is None (no minimum enforced).
        max_depth : int, optional
            The maximum depth of the phylogeny to prevent over-fitting in noisy data.
            Default is None (unlimited depth).
        min_clone_size : int, default=5
            The minimum number of cells required in a clone for it to be split further. Larger
            values help avoid spurious clones arising from noise in scRNA-seq-inferred CNAs.
            Recommended: 50-100 cells for noisy data, 10-50 for high-quality data.
        terminate_by : {'probability', 'BIC'}, default='probability'
            The criterion used to terminate clone splitting:
            - 'probability': Uses assignment confidence to handle uncertainty in noisy data
            - 'BIC': Uses Bayesian Information Criterion for model selection
        assignment_confidence_threshold : float, default=0.75
            Minimum confidence threshold for clone assignments when terminate_by='probability'.
            Higher values (0.8-0.9) recommended for very noisy scRNA-seq data to ensure
            confident assignments. Must be between 0 and 1.
        assignment_confidence_proportion : float, default=0.8
            Minimum proportion of cells with confident assignments required for clone splitting
            when terminate_by='probability'. Higher values help avoid splitting based on
            uncertain assignments in noisy data. Must be between 0 and 1.
        bic_penalty_strength : float, default=1.0
            Strength of BIC penalty term. Higher values (>1.0) encourage simpler models,
            useful for noisy data to prevent over-fitting.

        Attributes
        ----------
        terminal_clones : Dict[str, pd.Index]
            Dictionary tracking clones marked as terminal (no further splitting).
            Keys are clone identifiers, values are pandas Index objects of cell identifiers.
        clones : Dict[str, pd.Index]
            Dictionary mapping current clone IDs to pandas Index objects of
            cell identifiers belonging to each clone. Updated during tree construction.
        depth : int
            Current depth of phylogenetic tree construction.

        Raises
        ------
        AssertionError
            If character_matrix is not a pandas DataFrame.
        ValueError
            If character_matrix cannot be converted to integer values.
        AssertionError
            If confidence thresholds are not between 0 and 1, or if min/max depth values
            are invalid.

        Examples
        --------
        Basic usage with scRNA-seq-inferred CNA data:

        >>> from picasso import Picasso, load_data
        >>>
        >>> # Load example CNA data
        >>> character_matrix = load_data()
        >>>
        >>> # Initialize PICASSO with parameters suitable for noisy data
        >>> picasso = Picasso(character_matrix,
        ...                  min_clone_size=10,  # Choose a larger value for very noisy data
        ...                  assignment_confidence_threshold=0.85,  # Higher confidence
        ...                  assignment_confidence_proportion=0.9)
        >>>
        >>> # Fit the model
        >>> picasso.fit()
        >>>
        >>> # Get results
        >>> phylogeny = picasso.get_phylogeny()
        >>> clone_assignments = picasso.get_clone_assignments()

        For very noisy data, use stricter parameters:

        >>> # Parameters for very noisy scRNA-seq-inferred CNAs
        >>> picasso_strict = Picasso(character_matrix,
        ...                         min_clone_size=50,
        ...                         max_depth=8,  # Limit depth to avoid over-fitting
        ...                         assignment_confidence_threshold=0.9,
        ...                         assignment_confidence_proportion=0.95)  # Stronger penalty
        >>> # Alternatively, use BIC-based termination
        >>> picasso_strict = Picasso(character_matrix,
        ...                         min_clone_size=50,
        ...                         min_depth=3, # Force splitting to a depth of 3
        ...                         max_depth=8,  # Limit depth to avoid over-fitting
        ...                         terminate_by='BIC')
        >>> picasso_strict.fit()

        Notes
        -----
        The PICASSO algorithm proceeds through the following steps:

        1. **Initialization**: All cells start in a single root clone
        2. **Iterative Splitting**: At each depth level:
           - For each current clone, fit Categorical Mixture Models with k=1 and k=2 components
           - Evaluate splitting criteria (BIC or assignment confidence)
           - Split clones that meet criteria into two daughter clones
        3. **Termination**: Stop when no clones can be split further or max_depth is reached
        4. **Tree Construction**: Build phylogenetic tree from clone hierarchy. Leaves are clones containing cells
            whose CNAs cannot be further distinguised reliably.

        **Handling Noisy scRNA-seq Data**:
        - Uses probabilistic assignment with confidence thresholds
        - Minimum clone size prevents spurious small clones from noise
        - BIC penalty prevents over-fitting to noise artifacts
        - Confidence-based termination handles assignment uncertainty

        **Model Assumptions**:
        - CNAs are acquired progressively but can be acquired multiple times independently (no perfect phylogeny assumption)
        - Each genomic feature evolves independently
        - Copy number states follow categorical distributions within clones
        - Noise is handled through mixture model uncertainty quantification

        See Also
        --------
        CloneTree : Class for phylogenetic tree visualization and analysis
        get_phylogeny : Method to extract the reconstructed phylogeny
        get_clone_assignments : Method to get cell-to-clone assignments

        """
        assert isinstance(
            assignment_confidence_threshold, float
        ), "assignment_confidence_threshold must be a float"
        assert isinstance(
            assignment_confidence_proportion, float
        ), "assignment_confidence_proportion must be a float"
        assert (
            0 <= assignment_confidence_threshold <= 1
        ), "assignment_confidence_threshold must be between 0 and 1"
        assert (
            0 <= assignment_confidence_proportion <= 1
        ), "assignment_confidence_proportion must be between 0 and 1"
        self.assignment_confidence_threshold = assignment_confidence_threshold
        self.assignment_confidence_proportion = assignment_confidence_proportion
        self.bic_penalty_strength = bic_penalty_strength

        assert isinstance(
            character_matrix, pd.DataFrame
        ), "character_matrix must be a pandas DataFrame"
        # Convert character matrix to integer values
        try:
            character_matrix = character_matrix.astype(int).copy()
        except:
            raise ValueError("Character matrix must be convertible to integer values.")

        assert (
            isinstance(min_depth, int) or min_depth is None
        ), "min_depth must be an integer or None"
        assert (
            isinstance(max_depth, int) or max_depth is None
        ), "max_depth must be an integer or None"
        assert (
            isinstance(min_clone_size, int) or min_clone_size is None
        ), "min_clone_size must be an integer or None"
        terminate_by = terminate_by.upper()
        assert terminate_by in [
            "PROBABILITY",
            "BIC",
        ], 'terminate_by must be either "probability" or "BIC"'

        self.character_matrix = character_matrix
        self.min_depth = min_depth if min_depth is not None else 0
        self.max_depth = max_depth if max_depth is not None else float("inf")
        if min_clone_size is not None:
            assert isinstance(min_clone_size, int), "min_clone_size must be an integer"
            assert min_clone_size > 0, "min_clone_size must be greater than 0"
        self.min_clone_size = min_clone_size if min_clone_size is not None else 1
        self.terminate_by = terminate_by

        self.terminal_clones = {}
        self.clones = {"1": character_matrix.index}
        self.depth = 0

    def split_clone(self, clone: str, force_split: bool = False) -> Dict[str, pd.Index]:
        """
        Attempt to split a single clone into two daughter clones using mixture modeling.

        Evaluates whether a clone should be split by fitting Categorical Mixture Models
        and applying termination criteria. This is the core method for handling noisy
        CNA data through probabilistic modeling and confidence-based decisions.

        Parameters
        ----------
        clone : str
            Identifier of the clone to attempt splitting. Should be a key in self.clones.
        force_split : bool, default=False
            If True, override normal termination criteria and force splitting (used when
            min_depth hasn't been reached). Still respects minimum clone size constraints.

        Returns
        -------
        dict
            Dictionary mapping new clone identifiers to pandas Index objects containing
            the cell/sample identifiers assigned to each clone:
            - If split successful: {'{clone}-0': cells_0, '{clone}-1': cells_1}
            - If terminated: {'{clone}-STOP': original_cells}
            - If already terminal: {clone: original_cells}

        Examples
        --------
        >>> from picasso import Picasso, load_data
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix)
        >>> # After some fitting steps, try splitting a specific clone
        >>> result = picasso.split_clone('1-0')
        >>> print(f"Split result: {list(result.keys())}")

        Force splitting (ignoring confidence criteria):
        >>> forced_result = picasso.split_clone('1-1', force_split=True)

        Notes
        -----
        **Splitting Process**:
        1. Check if clone is already terminal (return unchanged)
        2. Extract CNA profiles for cells in the clone
        3. Filter features with sufficient variance (> 1e-10) for performance improvements
        4. Fit mixture models with k=1 and k=2 components
        5. Evaluate termination criteria (BIC or confidence)
        6. Apply minimum clone size constraint
        7. Return split result or mark as terminal

        **Termination Criteria**:
        - **BIC**: k=1 model has lower BIC than k=2 model
        - **Probability**: Insufficient assignment confidence or proportion
        - **Size constraint**: Either daughter clone below min_clone_size

        **Noise Handling**:
        - Confidence thresholds prevent splits based on uncertain assignments
        - Minimum clone sizes avoid spurious small clusters
        - Variance filtering removes uninformative features
        - Multiple model fitting attempts with different initializations

        See Also
        --------
        step : Apply split_clone to all current leaf clones
        _select_model : Internal method for mixture model fitting
        fit : Main method that orchestrates the complete splitting process
        """
        new_clones = {}
        log.debug(f"\t Processing Clone {clone} of size {len(self.clones[clone])}.")
        if clone in self.terminal_clones:
            new_clones[clone] = self.clones[clone]
            return new_clones

        # Get the samples corresponding to this leaf node
        samples = self.clones[clone]
        # Get the character matrix for these samples, keeping only features with variance greater than 1e-10
        character_matrix = self.character_matrix.loc[
            samples, self.character_matrix.var() > 1e-10
        ].copy()

        # Ensure that the character matrix is not empty
        if len(character_matrix.columns) == 0:
            # Clone is terminal and cannot be split further
            self.terminal_clones[clone] = samples
            new_clones[clone] = samples
            return new_clones

        # Ensure that the character matrix has integer values with a minimum of 0
        X = copy.deepcopy(character_matrix.values - character_matrix.min().min())
        # Fit a Categorical Mixture Model to the character matrix
        model2, bic2 = self._select_model(X, 2)
        # Split the clone into two children
        responsibilities = model2.predict_proba(X).numpy()
        assignments = np.argmax(responsibilities, axis=1)

        terminate = False
        if self.terminate_by == "BIC":
            model1, bic1 = self._select_model(X, 1)
            if bic1 < bic2:
                terminate = True

        if self.terminate_by == "PROBABILITY":
            # Determine confident assignments
            confident_assignments = (
                np.max(responsibilities, axis=1) >= self.assignment_confidence_threshold
            )
            confident_proportion = (
                np.sum(confident_assignments) / character_matrix.shape[0]
            )
            if confident_proportion < self.assignment_confidence_proportion:
                terminate = True

        if self.terminate_by == "CHI_SQUARED":
            terminate = not self._perform_chi_squared(
                X, responsibilities, self.chi_squared_p_value
            )

        try:
            samples_in_clone_0 = samples[assignments == 0]
            samples_in_clone_1 = samples[assignments == 1]
        except Exception as e:
            print(set(assignments))
            print(assignments == 0)
            raise e

        # If the algorithm is forced to split, try to split the clone regardless of the BIC score
        if force_split and terminate:
            log.debug(f"\t -Forced split of clone {clone}.")
            terminate = False

        # No matter what, if a clone is too small, terminate it
        if (
            len(samples_in_clone_0) < self.min_clone_size
            or len(samples_in_clone_1) < self.min_clone_size
        ):
            terminate = True

        if terminate:
            self.terminal_clones[f"{clone}-STOP"] = samples
            new_clones[f"{clone}-STOP"] = samples
            log.debug(f"\t -Terminated clone {clone}.")

        else:
            new_clones[f"{clone}-0"] = samples_in_clone_0
            new_clones[f"{clone}-1"] = samples_in_clone_1
            log.debug(
                f'\t -Split clone {clone} into sublones of sizes {len(new_clones[f"{clone}-0"])}'
                f' and {len(new_clones[f"{clone}-1"])}'
            )
        log.debug("\t -------------------")
        return new_clones

    def step(self, force_split: bool = False) -> None:
        """
        Execute one complete iteration of clone splitting across all current leaf clones.

        Applies the split_clone method to all current leaf clones in parallel, representing
        one depth level of the phylogenetic reconstruction process. This method coordinates
        the simultaneous evaluation of all clones at the current tree depth.

        Parameters
        ----------
        force_split : bool, default=False
            If True, attempts to force splits even when normal termination criteria
            are met. Used when enforcing minimum tree depth requirements. Individual
            clones may still be terminated if size constraints are violated.

        Notes
        -----
        **Single Step Process**:
        1. Iterate through all current leaf clones
        2. Apply split_clone() to each clone
        3. Collect all resulting clones (split or terminal)
        4. Update self.clones with the new clone structure
        5. Terminal clones are tracked in self.terminal_clones

        **Progress Tracking**:
        - Uses tqdm progress bar to show splitting progress
        - Logs clone processing information at debug level
        - Reports clone sizes and splitting decisions

        **State Modification**:
        - Updates self.clones with new clone structure
        - Adds terminal clones to self.terminal_clones
        - Preserves cell-to-clone assignment mappings

        **Parallelization Note**:
        Currently processes clones sequentially. Future versions may implement
        parallel processing for large datasets.

        Examples
        --------
        >>> from picasso import Picasso, load_data
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix)
        >>> print(f"Initial clones: {len(picasso.clones)}")
        >>> picasso.step()  # Perform one splitting iteration
        >>> print(f"After step: {len(picasso.clones)} clones, {len(picasso.terminal_clones)} terminal")

        Force splitting to explore deeper structure:
        >>> picasso.step(force_split=True)

        See Also
        --------
        split_clone : Method applied to individual clones during this step
        fit : Complete algorithm that calls step() iteratively until termination
        """
        new_clones = {}
        for clone in tqdm(self.clones):
            # Get the size of the clone
            updated_clones = self.split_clone(clone, force_split)
            for key, value in updated_clones.items():
                new_clones[key] = value
        self.clones = new_clones

    def fit(self) -> None:
        """
        Fit the PICASSO phylogenetic model to the noisy CNA data.

        Executes the complete PICASSO algorithm by iteratively splitting clones until
        termination criteria are met. The algorithm is designed to handle noise and
        uncertainty in scRNA-seq-inferred CNA data through probabilistic modeling and
        confidence-based termination.

        Parameters
        ----------
        None

        Returns
        -------
        None
            Modifies the instance in-place by updating clones, terminal_clones, and depth.

        Notes
        -----
        The fitting process proceeds as follows:

        1. **Iterative Splitting**: At each depth level, all current leaf clones are
           evaluated for splitting using Categorical Mixture Models
        2. **Noise Handling**: Uses confidence thresholds and minimum clone sizes to
           avoid splits driven by noise artifacts
        3. **Forced Splitting**: If min_depth is specified, forces splits until that
           depth is reached (unless clone size is insufficient)
        4. **Termination**: Stops when all clones are terminal, max_depth is reached,
           or no clones meet splitting criteria

        **Termination Conditions**:
        - All leaf clones have been marked as terminal
        - Maximum depth limit reached (if specified)
        - No clones have sufficient size for splitting
        - Confidence/BIC criteria not met for any clones

        **For Noisy scRNA-seq Data**:
        - Higher confidence thresholds prevent spurious splits
        - Larger minimum clone sizes reduce noise-driven artifacts
        - BIC penalty helps prevent over-fitting to noise

        Examples
        --------
        >>> from picasso import Picasso, load_data
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix, min_clone_size=8)
        >>> picasso.fit()  # Fit the model
        >>> print(f"Final tree depth: {picasso.depth}")
        >>> print(f"Number of terminal clones: {len(picasso.terminal_clones)}")

        See Also
        --------
        step : Perform a single splitting iteration
        split_clone : Split an individual clone
        get_phylogeny : Extract the fitted phylogenetic tree
        """
        algorithm_finished = False
        start_time = time.time()
        while not algorithm_finished:
            self.depth += 1
            log.info(
                f"Tree Depth {self.depth}: {len(self.clones)} clone(s), {len(self.terminal_clones)} terminal clone(s). "
                f"Force Split: {self.depth <= self.min_depth}"
            )
            self.step(force_split=self.depth <= self.min_depth)

            # Determine whether all leaf nodes have been terminated or if the algorithm has reached the maximum depth
            if self.depth < self.min_depth:
                continue
            if self.depth >= self.max_depth:
                log.info(f"Maximum depth of {self.max_depth} reached.")
                algorithm_finished = True
            if len(set(self.clones)) == len(set(self.terminal_clones)):
                log.info("All leaf nodes have been terminated.")
                algorithm_finished = True
        log.info(
            f"PICASSO algorithm finished in {time.time() - start_time:.2f} seconds."
        )

    def _select_model(
        self, X: np.ndarray, n_clusters: int
    ) -> Tuple[GeneralMixtureModel, float]:
        """
        Select the best Categorical Mixture Model using multiple random initializations.

        Fits multiple Categorical Mixture Models with different random initializations
        and selects the model with the best (lowest) BIC score. This approach helps
        overcome local optima that can occur when fitting mixture models to noisy
        scRNA-seq-inferred CNA data.

        Parameters
        ----------
        X : np.ndarray
            Character matrix for a clone with cells as rows and genomic features as columns.
            Values should be non-negative integers representing copy number states.
        n_clusters : int
            Number of mixture components (clusters) to fit. Must be positive.
            For n_clusters=1, fits a single Categorical distribution.
            For n_clusters>1, fits a mixture of Categorical distributions.

        Returns
        -------
        tuple
            Tuple containing:
            - best_model : GeneralMixtureModel
                Fitted mixture model with the lowest BIC score across all trials
            - best_bic : float
                BIC score of the best model (lower values indicate better fit)

        Notes
        -----
        **Multi-trial Strategy**:
        - For single components (n_clusters=1): Runs 1 trial (deterministic)
        - For multiple components: Runs up to 5 trials with different initializations
        - Selects model with lowest BIC score across all successful trials

        **Error Handling**:
        - If standard random initialization fails, tries controlled initialization
        - Uses _initialize_clusters() method for robust fallback initialization
        - Logs trial information and BIC scores for debugging

        **BIC Score Calculation**:
        - Uses custom _get_BIC_score() method with configurable penalty strength
        - Accounts for model complexity and data likelihood
        - Lower BIC scores indicate better balance of fit and simplicity

        Raises
        ------
        AssertionError
            If BIC score is infinite or NaN, indicating model fitting failure.

        See Also
        --------
        _get_BIC_score : Method for calculating BIC scores
        _initialize_clusters : Method for controlled model initialization
        split_clone : Method that uses this function for clone splitting decisions
        """
        n_trials = 0
        best_bic = np.inf
        best_model = None

        if n_clusters == 1:
            max_trials = 1
        else:
            max_trials = 5

        while n_trials < max_trials:
            try:
                if n_clusters == 1:
                    distributions = [Categorical().fit(X)]
                else:
                    distributions = [Categorical() for _ in range(n_clusters)]
                model = GeneralMixtureModel(distributions, verbose=False).fit(X)
                bic_score = self._get_BIC_score(model, X, self.bic_penalty_strength)
                assert not np.isinf(bic_score), f"BIC score is {bic_score}."
                if n_clusters > 1:
                    log.debug(f"\t -Trial {n_trials + 1}: BIC = {bic_score}")
                # Select the lowest BIC score
                if bic_score < best_bic:
                    best_bic = bic_score
                    best_model = model
                n_trials += 1
            except Exception as e:
                log.debug(
                    f"\t -Trial {n_trials + 1}: {e}. Retrying with controlled initialization."
                )
                distributions = self._initialize_clusters(X, n_clusters)
                model = GeneralMixtureModel(distributions, verbose=False).fit(X)
                bic_score = self._get_BIC_score(model, X)
                assert not np.isinf(bic_score), f"BIC score is {bic_score}."
                if n_clusters > 1:
                    log.debug(f"\t -Trial {n_trials + 1}: BIC = {bic_score}")
                # Select the lowest BIC score
                if bic_score < best_bic:
                    best_bic = bic_score
                    best_model = model
                n_trials += 1
        assert not np.isnan(best_bic) and not np.isinf(
            best_bic
        ), f"Best BIC score is {best_bic}."
        return best_model, best_bic

    @staticmethod
    def _get_BIC_score(
        model: GeneralMixtureModel, X: np.ndarray, bic_penalty_strength: float = 1.0
    ) -> float:
        """
        Calculate Bayesian Information Criterion (BIC) score for mixture model selection.

        Computes the BIC score to balance model fit quality against complexity. Lower BIC
        scores indicate better models. The BIC penalizes complex models to prevent over-fitting
        to noise in scRNA-seq-inferred CNA data.

        Parameters
        ----------
        model : GeneralMixtureModel
            Fitted Categorical Mixture Model from pomegranate library.
        X : np.ndarray
            Input data matrix used for model fitting, with cells as rows and
            genomic features as columns.
        bic_penalty_strength : float, default=1.0
            Strength of the complexity penalty term. Values >1.0 encourage simpler
            models, which can help prevent over-fitting in noisy data.

        Returns
        -------
        float
            BIC score where lower values indicate better models. Calculated as:
            BIC = -2 * log_likelihood + penalty_strength * n_params * log(n_samples)

        Notes
        -----
        **BIC Formula Components**:
        - Log-likelihood: Measures how well the model fits the data
        - Complexity penalty: Number of model parameters × log(sample size)
        - Penalty strength: Allows adjustment of complexity penalty weight

        **Parameter Counting**:
        - For each cluster: (n_states - 1) × n_features parameters
        - Plus (n_clusters - 1) mixture weight parameters
        - Total complexity scales with number of clusters and features

        **Interpretation**:
        - More negative log probability → higher BIC (worse)
        - More model parameters → higher BIC (worse)
        - Optimal model minimizes BIC score

        **Noise Handling**:
        - Higher penalty_strength values help avoid over-fitting to noise
        - Particularly important for noisy scRNA-seq-inferred CNAs

        See Also
        --------
        _select_model : Method that uses BIC scores for model selection
        split_clone : Method that compares BIC scores to decide on splitting
        """
        D, K = model.distributions[0].probs.shape
        params_per_cluster = (
            model.distributions[0].probs.shape[1] - 1
        ) * model.distributions[0].probs.shape[0]
        n_clusters = len(model.distributions)
        n_params = params_per_cluster * n_clusters + n_clusters - 1
        logprob = model.log_probability(X).sum()
        bic_score = -2 * logprob + bic_penalty_strength * (
            n_params * np.log(X.shape[0])
        )

        return bic_score

    @staticmethod
    def _perform_chi_squared(
        X: np.ndarray, responsibilities: np.ndarray, threshold: float = 0.05
    ) -> bool:
        """
        Perform chi-squared test to evaluate statistical significance of clone splitting.

        Tests whether the observed differences in CNA profiles between two potential
        clones are statistically significant using a chi-squared test of independence.
        This provides an alternative to BIC-based splitting criteria.

        Parameters
        ----------
        X : np.ndarray
            Character matrix with cells as rows and genomic features as columns.
            Values should be non-negative integers representing copy number states.
        responsibilities : np.ndarray
            Soft assignment probabilities from mixture model with shape (n_cells, 2).
            Each row sums to 1.0 and represents the probability of each cell belonging
            to each of the two potential clones.
        threshold : float, default=0.05
            Significance threshold for the chi-squared test. Lower values require
            stronger evidence for splitting.

        Returns
        -------
        bool
            True if the chi-squared test supports splitting (p < threshold).
            False if splitting is not statistically supported.

        Notes
        -----
        **Statistical Test**:
        - Constructs contingency tables of copy number states × clones
        - Uses weighted counts based on soft assignment probabilities
        - Performs chi-squared test of independence

        **Interpretation**:
        - Low p-value (< threshold): Significant difference between clones → split
        - High p-value (≥ threshold): No significant difference → don't split

        **Implementation Details**:
        - Handles all unique copy number states in the data
        - Uses soft assignments to weight expected frequencies
        - Adds small pseudocount (1e-5) to avoid zero frequencies
        - Flattens contingency tables for chi-squared test

        **Usage Context**:
        - Alternative to BIC-based termination criteria
        - Can be selected via terminate_by='CHI_SQUARED' parameter
        - Provides statistical foundation for splitting decisions

        See Also
        --------
        split_clone : Method that can use this test for termination decisions
        scipy.stats.chi2_contingency : Underlying statistical test function
        """
        from scipy.stats import chi2_contingency, chisquare

        # States are assumed to be 0, 1, 2, ... (already shifted to be positive)
        states = np.unique(X)

        # Number of regions
        num_regions = X.shape[1]

        # Initialize the expected frequencies table
        expected_frequencies = np.zeros((len(states), num_regions, 2))

        # Calculate weighted frequencies for each (state, region) pair in each clone using numpy broadcasting
        for clone in range(responsibilities.shape[1]):
            for state in states:
                # Get the cells with the current state
                state_cells = X == state
                # Get the expected frequency for the current state in the current clone
                expected_frequencies[state, :, clone] = (
                    state_cells.T @ responsibilities[:, clone]
                )

        # Ensure that there are no expected frequencies of zero
        expected_frequencies[expected_frequencies == 0] = 1e-5

        # For performing the Chi-squared test, flatten the tables appropriately
        contingency_table_clone1 = expected_frequencies[:, :, 0].flatten()
        contingency_table_clone2 = expected_frequencies[:, :, 1].flatten()
        contingency_table = np.vstack(
            [contingency_table_clone1, contingency_table_clone2]
        )

        # Perform Chi-squared test
        chi2, p, dof, expected = chi2_contingency(contingency_table)

        # Decision based on the p-value
        log.debug("Chi-squared p-value:", p)
        split_decision = p < threshold
        return split_decision

    @staticmethod
    def _initialize_clusters(
        array: np.ndarray, num_partitions: int
    ) -> List[Categorical]:
        """
        Initialize Categorical distributions for mixture model with controlled partitioning.

        Creates robust initial parameter estimates for Categorical mixture models by
        randomly partitioning cells and fitting separate distributions to each partition.
        This controlled initialization helps overcome convergence issues that can occur
        with random initialization on noisy CNA data.

        Parameters
        ----------
        array : np.ndarray
            Character matrix with cells as rows and genomic features as columns.
            Values should be non-negative integers representing copy number states.
        num_partitions : int
            Number of mixture components to initialize. Must be at least 2.

        Returns
        -------
        List[Categorical]
            List of fitted Categorical distributions, one for each partition.
            Each distribution is pre-fitted to a subset of the input data.

        Raises
        ------
        ValueError
            If num_partitions is less than 2.

        Notes
        -----
        **Initialization Strategy**:
        - Randomly permutes cell indices to avoid bias
        - Creates random breakpoints for partitioning
        - Ensures each partition has at least one cell
        - Augments each partition with maximum values to ensure full state coverage

        **Robustness Features**:
        - Adds row with maximum values to each partition
        - Ensures all copy number states are represented in each distribution
        - Prevents model fitting failures due to missing states
        - Creates diverse initial parameter estimates

        **Use Case**:
        - Fallback initialization when standard random initialization fails
        - Called by _select_model() when mixture model fitting encounters errors
        - Particularly useful for noisy or sparse CNA data

        **Algorithm Steps**:
        1. Shuffle cell indices randomly
        2. Choose random breakpoints for partitioning
        3. Split indices into partitions
        4. Map indices back to original data
        5. Augment each partition with maximum values
        6. Fit Categorical distribution to each partition

        See Also
        --------
        _select_model : Method that uses this function as fallback initialization
        pomegranate.distributions.Categorical : Distribution class being initialized
        """
        if num_partitions < 2:
            raise ValueError("num_partitions must be at least 2")

        # Shuffle the array indices
        num_elements = array.shape[0]
        shuffled_indices = np.random.permutation(num_elements)

        # Choose random break points
        break_points = sorted(
            np.random.choice(num_elements - 1, num_partitions - 1, replace=False) + 1
        )

        # Split the indices at the break points
        partitions = np.split(shuffled_indices, break_points)

        # Map indices back to the original array
        partitioned_arrays = [array[partition] for partition in partitions]

        # For each partition, add one row containing the maximum value of the overall array in each column
        max_values = np.max(array, axis=0)
        for i in range(num_partitions):
            partitioned_arrays[i] = np.vstack([partitioned_arrays[i], max_values])

        distributions = [Categorical().fit(y) for y in partitioned_arrays]

        return distributions

    def get_phylogeny(self) -> ete3.Tree:
        """
        Extract the reconstructed phylogenetic tree from the fitted PICASSO model.

        Converts the hierarchical clone structure into an ete3.Tree object for
        visualization and downstream analysis. The tree represents the inferred
        evolutionary relationships between clones based on their CNA profiles.

        Returns
        -------
        ete3.Tree
            Phylogenetic tree where leaves represent terminal clones and internal
            nodes represent ancestral clones. Node names correspond to clone IDs
            from the splitting process (e.g., '1', '1-0', '1-1', '1-0-STOP').

        Examples
        --------
        >>> from picasso import Picasso, load_data
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix)
        >>> picasso.fit()
        >>> tree = picasso.get_phylogeny()
        >>> print(tree.get_ascii())  # Display tree structure
        >>> print(f"Tree has {len(tree.get_leaves())} terminal clones")

        Get leaf names:
        >>> leaf_names = tree.get_leaf_names()
        >>> print(f"Terminal clones: {leaf_names}")

        Notes
        -----
        - The tree topology reflects the binary splitting process used by PICASSO
        - Internal nodes represent decision points where clones were split
        - Terminal nodes (leaves) represent final clones that could not be split further
        - Node names encode the splitting history (e.g., '1-0-1' = root -> left -> right)
        - Trees from noisy data may have different topologies due to uncertainty handling

        See Also
        --------
        get_clone_assignments : Get cell-to-clone assignments
        CloneTree : Class for enhanced tree visualization and analysis
        create_tree_from_paths : Static method for tree construction from paths
        """
        phylogeny = self.create_tree_from_paths(self.clones.keys(), "-")
        return phylogeny

    def get_clone_assignments(self) -> pd.DataFrame:
        """
        Extract cell-to-clone assignments from the fitted PICASSO model.

        Returns a DataFrame mapping each cell/sample to its assigned terminal clone.
        These assignments represent the final clustering result after the phylogenetic
        reconstruction process.

        Returns
        -------
        pd.DataFrame
            DataFrame with cell/sample identifiers as index and a 'clone_id' column
            containing the assigned clone ID for each cell. Clone IDs correspond to
            the terminal nodes in the phylogenetic tree.

        Examples
        --------
        >>> from picasso import Picasso, load_data
        >>> character_matrix = load_data()
        >>> picasso = Picasso(character_matrix)
        >>> picasso.fit()
        >>> assignments = picasso.get_clone_assignments()
        >>> print(assignments.head())
        >>> print(f"Number of clones: {assignments['clone_id'].nunique()}")

        Get cells in a specific clone:
        >>> clone_cells = assignments[assignments['clone_id'] == '1-0-STOP'].index
        >>> print(f"Cells in clone 1-0-STOP: {list(clone_cells)}")

        Clone size distribution:
        >>> clone_sizes = assignments['clone_id'].value_counts()
        >>> print("Clone sizes:")
        >>> print(clone_sizes)

        Notes
        -----
        - Each cell is assigned to exactly one terminal clone
        - Clone IDs reflect the splitting hierarchy (e.g., '1-0-STOP', '1-1-0-STOP')
        - The '-STOP' suffix indicates terminal clones that were not split further
        - Assignment quality depends on the noise level in the input CNA data
        - For very noisy data, some assignments may have lower confidence

        See Also
        --------
        get_phylogeny : Get the phylogenetic tree structure
        CloneTree : Class for integrated analysis of assignments and phylogeny
        fit : Method that performs the clustering and phylogeny reconstruction
        """
        clone_assigments = {"samples": [], "clone_id": []}
        for clone in self.clones:
            clone_assigments["samples"].extend(self.clones[clone])
            clone_assigments["clone_id"].extend([clone] * len(self.clones[clone]))
        clone_assigments = pd.DataFrame(clone_assigments).set_index("samples")
        return clone_assigments

    @staticmethod
    def create_tree_from_paths(paths: List[str], separator: str = ":") -> ete3.TreeNode:
        """
        Construct phylogenetic tree from hierarchical clone path identifiers.

        Converts a list of clone path strings into an ete3 tree structure by parsing
        the hierarchical splitting history encoded in each path. This is used internally
        by PICASSO to generate the final phylogenetic tree from the clone splitting process.

        Parameters
        ----------
        paths : list of str
            List of clone path identifiers representing the hierarchical structure.
            Each path encodes the splitting history (e.g., '1', '1-0', '1-0-STOP', '1-1-0').
            All paths must start with the same root character.
        separator : str, default=':'
            Character used to separate levels in the path hierarchy. PICASSO uses '-'
            by default for clone paths.

        Returns
        -------
        ete3.TreeNode
            Root node of the constructed phylogenetic tree where:
            - Leaves represent terminal clones
            - Internal nodes represent ancestral states/splitting points
            - Node names correspond to the original path identifiers

        Examples
        --------
        Basic tree construction from clone paths:

        >>> from picasso.build_tree import Picasso
        >>>
        >>> # Example clone paths from PICASSO splitting
        >>> clone_paths = ['1', '1-0-STOP', '1-1-0-STOP', '1-1-1-STOP']
        >>> tree = Picasso.create_tree_from_paths(clone_paths, '-')
        >>> print(tree.get_ascii())
        >>> print(f"Leaves: {tree.get_leaf_names()}")

        Notes
        -----
        **Path Structure**:
        - Root level: Single character (typically '1')
        - Subsequent levels: Added via separator (e.g., '1-0', '1-1')
        - Terminal indicator: Often ends with '-STOP' for final clones

        **Tree Construction Logic**:
        - Identifies common root from all paths
        - Builds tree level by level based on path prefixes
        - Creates parent-child relationships following path hierarchy
        - Handles variable depth paths automatically

        **Internal Use**:
        - Called by get_phylogeny() to convert clone structure to tree
        - Maintains clone ID information in node names
        - Preserves splitting history for downstream analysis

        Raises
        ------
        AssertionError
            If paths don't share a common root character.

        See Also
        --------
        get_phylogeny : Public method that uses this function to create phylogenetic trees
        fit : Method that generates the clone paths through iterative splitting
        """
        paths = list(set(paths))
        root = list(set([path[0] for path in paths]))
        assert len(root) == 1, "All paths must start with the same character"
        max_depth = max([len(path.split(separator)) for path in paths])

        all_nodes = {str(root[0]): ete3.TreeNode(name=str(root[0]))}
        for depth in range(2, max_depth + 1):
            prefix_paths = []
            for path in paths:
                if len(path.split(separator)) < depth:
                    continue
                prefix_paths.append(separator.join(path.split(separator)[:depth]))
            prefix_paths = set(prefix_paths)
            for path in prefix_paths:
                parent = all_nodes[separator.join(path.split(separator)[:-1])]
                node = ete3.TreeNode(name=path)
                parent.add_child(node)
                all_nodes[path] = node

        return all_nodes[str(root[0])]


# Define public API
__all__ = ["Picasso"]
