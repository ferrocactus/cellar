import unittest

import numpy as np
from numpy.testing import (assert_almost_equal, assert_array_almost_equal,
                           assert_array_equal, assert_equal)
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from src.units._markers import _ttest_differential_expression


class ttest(unittest.TestCase):
    def test_ttest_differential_expression(self):
        alpha = 0.05
        markers_n = 5
        correction = 'hs'

        x1 = np.random.random(size=(7, 20))
        test_results1 = _ttest_differential_expression(
            x1, x1, alpha=alpha, markers_n=markers_n, correction=correction)
        assert_array_equal(np.array([]), test_results1['indices'])

        x21 = np.ones(shape=(4, 20)) @ np.diag(np.arange(27, 7, -1))
        x22 = np.ones(shape=(3, 20)) * (-4)
        test_results2 = _ttest_differential_expression(
            x21, x22, alpha=alpha, markers_n=markers_n, correction=correction)
        assert_array_equal(np.arange(markers_n), test_results2['indices'])
        test_results3 = _ttest_differential_expression(
            x22, x21, alpha=alpha, markers_n=markers_n, correction=correction)
        assert_array_equal(np.array([]), test_results3['indices'])


        x31 = np.array(
            [[-1, -0.5, 1, 5],
             [-0.5, -1, 0.5, 5],
             [0, 0, 1, 7]]
        )
        x32 = np.array(
            [[-0.5, 13, 0.5, -2],
             [-2, 12, 0, -1],
             [1, 11, -1, -1]]
        )
        test_results4 = _ttest_differential_expression(
            x31, x32, alpha=alpha, markers_n=markers_n, correction=correction)
        assert_array_equal(np.array([3]), test_results4['indices'])
