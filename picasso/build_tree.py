from pomegranate.gmm import GeneralMixtureModel
from pomegranate.distributions import Categorical

import copy
import numpy as np
import pandas as pd

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
    def __init__(self,
                 character_matrix,
                 min_depth=None,
                 max_depth=None,
                 min_clone_size=5,
                 terminate_by='probability',
                 assignment_confidence_threshold=0.75,
                 assignment_confidence_proportion=0.8,
                 bic_penalty_strength=1.0,):
        """
        Initialize the PICASSO model.
        :param character_matrix: (pd.DataFrame) An integer matrix where rows are samples and columns are features.
        :param min_depth: (int) The minimum depth of the phylogeny.
        :param max_depth: (int) The maximum depth of the phylogeny.
        :param min_clone_size: (int) The minimum number of samples in a clone.
        :param terminate_by: (str) The criterion to use to terminate the algorithm. Either 'probability' or 'BIC'.
        :param assignment_confidence_threshold: (float) The minimum confidence threshold for clone assignments if
                                                terminate_by is 'probability'. Must be between 0 and 1.
        :param assignment_confidence_proportion: (float) The minimum proportion of samples with confident assignments
                                                 for a clone to be split if terminate_by is 'probability'. Must be
                                                 between 0 and 1.
        """
        assert isinstance(assignment_confidence_threshold, float), 'assignment_confidence_threshold must be a float'
        assert isinstance(assignment_confidence_proportion, float), 'assignment_confidence_proportion must be a float'
        assert 0 <= assignment_confidence_threshold <= 1, 'assignment_confidence_threshold must be between 0 and 1'
        assert 0 <= assignment_confidence_proportion <= 1, 'assignment_confidence_proportion must be between 0 and 1'
        self.assignment_confidence_threshold = assignment_confidence_threshold
        self.assignment_confidence_proportion = assignment_confidence_proportion
        self.bic_penalty_strength = bic_penalty_strength

        assert isinstance(character_matrix, pd.DataFrame), 'character_matrix must be a pandas DataFrame'
        # Convert character matrix to integer values
        try:
            character_matrix = character_matrix.astype(int).copy()
        except:
            raise ValueError("Character matrix must be convertible to integer values.")

        assert isinstance(min_depth, int) or min_depth is None, 'min_depth must be an integer or None'
        assert isinstance(max_depth, int) or max_depth is None, 'max_depth must be an integer or None'
        assert isinstance(min_clone_size, int) or min_clone_size is None, 'min_clone_size must be an integer or None'
        terminate_by = terminate_by.upper()
        assert terminate_by in ['PROBABILITY', 'BIC'], 'terminate_by must be either "probability" or "BIC"'

        self.character_matrix = character_matrix
        self.min_depth = min_depth if min_depth is not None else 0
        self.max_depth = max_depth if max_depth is not None else float('inf')
        if min_clone_size is not None:
            assert isinstance(min_clone_size, int), 'min_clone_size must be an integer'
            assert min_clone_size > 0, 'min_clone_size must be greater than 0'
        self.min_clone_size = min_clone_size if min_clone_size is not None else 1
        self.terminate_by = terminate_by

        self.terminal_clones = {}
        self.clones = {'1': character_matrix.index}
        self.depth = 0

    def split_clone(self, clone, force_split=False):
        """
        Split a clone into two subclones using a Categorical Mixture Model.
        :param clone: (str) The clone to split.
        :param force_split: (bool) If True, force the algorithm to split all leaf nodes.
        :return: new_clones: (dict) A dictionary of the new clones (mapping clone names to sample lists)
        :modifies: self.terminal_clones
        """
        new_clones = {}
        log.debug(f'\t Processing Clone {clone} of size {len(self.clones[clone])}.')
        if clone in self.terminal_clones:
            new_clones[clone] = self.clones[clone]
            return new_clones

        # Get the samples corresponding to this leaf node
        samples = self.clones[clone]
        # Get the character matrix for these samples, keeping only features with variance greater than 1e-10
        character_matrix = self.character_matrix.loc[samples, self.character_matrix.var() > 1e-10].copy()

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
        if self.terminate_by == 'BIC':
            model1, bic1 = self._select_model(X, 1)
            if bic1 < bic2:
                terminate = True

        if self.terminate_by == 'PROBABILITY':
            # Determine confident assignments
            confident_assignments = np.max(responsibilities, axis=1) >= self.assignment_confidence_threshold
            confident_proportion = np.sum(confident_assignments) / character_matrix.shape[0]
            if confident_proportion < self.assignment_confidence_proportion:
                terminate = True

        if self.terminate_by == 'CHI_SQUARED':
            terminate = not self._perform_chi_squared(X, responsibilities, self.chi_squared_p_value)

        try:
            samples_in_clone_0 = samples[assignments == 0]
            samples_in_clone_1 = samples[assignments == 1]
        except Exception as e:
            print(set(assignments))
            print(assignments==0)
            raise e

        # If the algorithm is forced to split, try to split the clone regardless of the BIC score
        if force_split and terminate:
            log.debug(f'\t -Forced split of clone {clone}.')
            terminate = False

        # No matter what, if a clone is too small, terminate it
        if len(samples_in_clone_0) < self.min_clone_size or len(samples_in_clone_1) < self.min_clone_size:
            terminate = True

        if terminate:
            self.terminal_clones[f'{clone}-STOP'] = samples
            new_clones[f'{clone}-STOP'] = samples
            log.debug(f'\t -Terminated clone {clone}.')

        else:
            new_clones[f'{clone}-0'] = samples_in_clone_0
            new_clones[f'{clone}-1'] = samples_in_clone_1
            log.debug(
                f'\t -Split clone {clone} into sublones of sizes {len(new_clones[f"{clone}-0"])}'
                f' and {len(new_clones[f"{clone}-1"])}')
        log.debug('\t -------------------')
        return new_clones

    def step(self, force_split=False):
        """
        Perform one splitting step of the PICASSO algorithm. Each non-terminal leaf node is split into two children.
        :param force_split: (bool) If True, force the algorithm to split all leaf nodes.
        :return: None.
        :modifies: self.clones, self.terminal_clones
        """
        new_clones = {}
        for clone in tqdm(self.clones):
            # Get the size of the clone
            updated_clones = self.split_clone(clone, force_split)
            for key, value in updated_clones.items():
                new_clones[key] = value
        self.clones = new_clones

    def fit(self):
        """
        Fit the PICASSO model by iteratively splitting leaf nodes until the algorithm terminates.
        Wrapper function for the step() method then performs logic to determine whether algorithm has terminated.
        If a minimum depth is specified, the algorithm will continue until the minimum depth is reached forcing splits
        UNLESS a clone is of insufficient size as determined by the min_clone_size parameter.
        :return: None.
        :modifies: self.clones, self.terminal_clones, self.depth
        """
        algorithm_finished = False
        start_time = time.time()
        while not algorithm_finished:
            self.depth += 1
            log.info(
                f'Tree Depth {self.depth}: {len(self.clones)} clone(s), {len(self.terminal_clones)} terminal clone(s). '
                f'Force Split: {self.depth <= self.min_depth}')
            self.step(force_split=self.depth<=self.min_depth)

            # Determine whether all leaf nodes have been terminated or if the algorithm has reached the maximum depth
            if self.depth < self.min_depth:
                continue
            if self.depth >= self.max_depth:
                log.info(f'Maximum depth of {self.max_depth} reached.')
                algorithm_finished = True
            if len(set(self.clones)) == len(set(self.terminal_clones)):
                log.info('All leaf nodes have been terminated.')
                algorithm_finished = True
        log.info(f'PICASSO algorithm finished in {time.time() - start_time:.2f} seconds.')

    def _select_model(self, X, n_clusters):
        """
        Select the best Categorical Mixture Model for a given number of clusters using the BIC score.
        Runs up to 5 trials to fit the model, retrying with controlled initialization if an error occurs.

        :param X: (np.ndarray) The character matrix for a clone.
        :param n_clusters: (int) The number of clusters to fit.
        :return: best_model: (GeneralMixtureModel) The best model.
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
                assert not np.isinf(bic_score), f'BIC score is {bic_score}.'
                if n_clusters > 1:
                    log.debug(f'\t -Trial {n_trials + 1}: BIC = {bic_score}')
                # Select the lowest BIC score
                if bic_score < best_bic:
                    best_bic = bic_score
                    best_model = model
                n_trials += 1
            except Exception as e:
                log.debug(f'\t -Trial {n_trials + 1}: {e}. Retrying with controlled initialization.')
                distributions = self._initialize_clusters(X, n_clusters)
                model = GeneralMixtureModel(distributions, verbose=False).fit(X)
                bic_score = self._get_BIC_score(model, X)
                assert not np.isinf(bic_score), f'BIC score is {bic_score}.'
                if n_clusters > 1:
                    log.debug(f'\t -Trial {n_trials + 1}: BIC = {bic_score}')
                # Select the lowest BIC score
                if bic_score < best_bic:
                    best_bic = bic_score
                    best_model = model
                n_trials += 1
        assert not np.isnan(best_bic) and not np.isinf(best_bic), f'Best BIC score is {best_bic}.'
        return best_model, best_bic

    @staticmethod
    def _get_BIC_score(model, X, bic_penalty_strength=1.0):
        """
        Compute BIC score for the model.  BIC score should be as low as possible. A more negative log prob translates
        roughly to a higher BIC score, and vice versa. We ideally want the most positive log prob we can get.
        :param model: (GeneralMixtureModel) Fitted model
        :param alpha: (float) alpha parameter for Dirichlet prior; less than 1 to encourage sparsity
        :return: (float) BIC score
        """
        D, K = model.distributions[0].probs.shape
        params_per_cluster = (model.distributions[0].probs.shape[1] - 1) * model.distributions[0].probs.shape[0]
        n_clusters = len(model.distributions)
        n_params = params_per_cluster * n_clusters + n_clusters - 1
        logprob = model.log_probability(X).sum()
        bic_score = -2 * logprob + bic_penalty_strength*(n_params * np.log(X.shape[0]))

        return bic_score

    @staticmethod
    def _perform_chi_squared(X, responsibilities, threshold = 0.05):
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
                expected_frequencies[state, :, clone] = (state_cells.T @ responsibilities[:, clone])

        # Ensure that there are no expected frequencies of zero
        expected_frequencies[expected_frequencies == 0] = 1e-5

        # For performing the Chi-squared test, flatten the tables appropriately
        contingency_table_clone1 = expected_frequencies[:, :, 0].flatten()
        contingency_table_clone2 = expected_frequencies[:, :, 1].flatten()
        contingency_table = np.vstack([contingency_table_clone1, contingency_table_clone2])

        # Perform Chi-squared test
        chi2, p, dof, expected = chi2_contingency(contingency_table)

        # Decision based on the p-value
        log.debug('Chi-squared p-value:', p)
        split_decision = p < threshold
        return split_decision

    @staticmethod
    def _initialize_clusters(array, num_partitions):
        """
        Initialize Categorical distribution clusters by randomly partitioning an array into a specified number of
        distinct arrays with random sizes. To avoid errors in fitting the model, each partition will contain at least
        one element and has values ranging from 0 to the maximum value in the array.

        :param array (np.ndarray): The input array to be partitioned.
        :param num_partitions (int): The number of partitions to create.
        :return List[np.ndarray]: A list of partitioned arrays.
        """
        if num_partitions < 2:
            raise ValueError("num_partitions must be at least 2")

        # Shuffle the array indices
        num_elements = array.shape[0]
        shuffled_indices = np.random.permutation(num_elements)

        # Choose random break points
        break_points = sorted(np.random.choice(num_elements - 1, num_partitions - 1, replace=False) + 1)

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

    def get_phylogeny(self):
        """
        Get the phylogeny of the clones as an ete3.Tree object.
        :return: phylogeny: (ete3.Tree) The phylogeny of the clones.
        """
        phylogeny = self.create_tree_from_paths(self.clones.keys(), '-')
        return phylogeny

    def get_clone_assignments(self):
        """
        Get the clone assignments for each sample as a DataFrame.
        :return: clone_assignments: (pd.DataFrame) A DataFrame where rows are samples and columns are clone IDs.
        """
        clone_assigments = {'samples':[], 'clone_id':[]}
        for clone in self.clones:
            clone_assigments['samples'].extend(self.clones[clone])
            clone_assigments['clone_id'].extend([clone] * len(self.clones[clone]))
        clone_assigments = pd.DataFrame(clone_assigments).set_index('samples')
        return clone_assigments

    @staticmethod
    def create_tree_from_paths(paths, separator=':'):
        """
        Create a tree from a list of paths representing leaves as sequences of characters separated by a separator.
        The root of the tree is the first character in each path and successive binary splits are made based on the
        characters in the path (0 for left child, 1 for right child).
        :param paths:
        :param separator:
        :return:
        """
        paths = list(set(paths))
        root = list(set([path[0] for path in paths]))
        assert len(root) == 1, 'All paths must start with the same character'
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