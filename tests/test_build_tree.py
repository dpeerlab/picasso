import unittest
import pandas as pd
import numpy as np

from picasso import Picasso

class TestPicasso(unittest.TestCase):

    def setUp(self):
        # Set up a simple character matrix for testing
        data = np.array([[0, 1, 2], [1, 0, 2], [2, 1, 0], [0, 2, 1]])
        self.character_matrix = pd.DataFrame(data, columns=['A', 'B', 'C'])
        self.picasso = Picasso(self.character_matrix, min_depth=1, max_depth=5, min_clone_size=1, terminate_by='BIC')

    def test_initialization(self):
        self.assertEqual(self.picasso.min_depth, 1)
        self.assertEqual(self.picasso.max_depth, 5)
        self.assertEqual(self.picasso.min_clone_size, 1)
        self.assertEqual(self.picasso.terminate_by, 'BIC')

    def test_invalid_initialization(self):
        with self.assertRaises(AssertionError):
            Picasso('not a dataframe')

    def test_split_clone(self):
        new_clones = self.picasso.split_clone('1')
        self.assertTrue(len(new_clones) > 0)

    def test_split_clone_force_split(self):
        new_clones = self.picasso.split_clone('1', force_split=True)
        self.assertTrue(len(new_clones) > 0)

    def test_step(self):
        self.picasso.step()
        self.assertTrue(len(self.picasso.clones) > 0)

    def test_fit(self):
        self.picasso.fit()
        self.assertTrue(self.picasso.depth > 0)
        self.assertTrue(len(self.picasso.terminal_clones) > 0)

    def test_select_model(self):
        X = np.array([[0, 1], [1, 0], [1, 1], [0, 0]])
        model, bic = self.picasso._select_model(X, 2)
        self.assertIsNotNone(model)
        self.assertTrue(bic > 0)

    def test_get_BIC_score(self):
        X = np.array([[0, 1], [1, 0], [1, 1], [0, 0]])
        model, _ = self.picasso._select_model(X, 2)
        bic = self.picasso._get_BIC_score(model, X)
        self.assertTrue(bic > 0)

    def test_initialize_clusters(self):
        array = np.array([[0, 1], [1, 0], [1, 1], [0, 0]])
        distributions = self.picasso._initialize_clusters(array, 2)
        self.assertTrue(len(distributions) == 2)

    def test_get_phylogeny(self):
        phylogeny = self.picasso.get_phylogeny()
        self.assertIsNotNone(phylogeny)

    def test_get_clone_assignments(self):
        clone_assignments = self.picasso.get_clone_assignments()
        self.assertIsInstance(clone_assignments, pd.DataFrame)
        self.assertTrue('clone_id' in clone_assignments.columns)

if __name__ == '__main__':
    unittest.main()
