import warnings

import numpy as np
import scipy as sp
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def welch_ttest(x1, x2):
    """
    Welch t-test for multiple arrays.

    Parameters
    __________
    x1: array, shape (n1, n_features)
    x2: array, shape (n2, n_features)

    Returns
    _______
    p: array, shape (n_features,)
        p-values for every column in x

    """
    eps = 1e-9
    n1 = x1[:, 0].size
    n2 = x2[:, 0].size
    m1 = np.mean(x1, axis=0)
    m2 = np.mean(x2, axis=0)
    v1 = np.var(x1, axis=0, ddof=1) + eps # add eps to avoid division by 0
    v2 = np.var(x2, axis=0, ddof=1) + eps # add eps to avoid division by 0
    t = (m1 - m2) / np.sqrt(v1 / n1 + v2 / n2)
    df = (v1 / n1 + v2 / n2)**2 / (v1**2 / (n1**2 * (n1 - 1))\
                + v2**2 / (n2**2 * (n2 - 1)))
    p = 2 * sp.stats.t.cdf(-abs(t), df)
    return p


def _ttest_differential_expression(x_i, x_not_i, alpha=0.05, max_n_genes=200,
                                   correction='hs'):
    """
    Find significant genes in pop1 compared to pop2 using ttest.

    Parameters
    __________

    x_i: array, shape (n, n_features)
        Target population. (Assumed not empty)

    x_not_i: array, shape (n_samples - n, n_features)
        The rest of the population. (Assumed not empty)

    alpha: float, between 0 and 1
        Alpha value to use for hypothesis testing.

    max_n_genes: int
        Number of top markers to return. Note: return number could be smaller.

    correction: string
        Correction method to apply for multitest.
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    _______

    test_results: dict
        Dictionary of the form {
            'indices': array, shape (<=max_n_genes,)
            'pvals': array, shape (<=max_n_genes,)
            'diffs': array, shape (<=max_n_genes,)
        }
        where indices correspond to the indices of the columns which
        were found to be significant, pvals corresponds to the p-value
        as returned by multipletests(ttest), and diffs corresponds to
        the different in the means of the vectors x_in vs x_not_in
        for the given significant columns.

    """
    if x_i.shape[0] <= 1 or x_not_i.shape[0] <= 1:
        warnings.warn(
            "Array contains <=1 elements, ttest results may be inaccurate.")

    pvals = welch_ttest(x_i, x_not_i)

    # difference of means of x_ij and x_not_ij
    diffs = np.mean(x_i, axis=0) - np.mean(x_not_i, axis=0)
    diffs[diffs <= 0] = -np.Inf

    if correction != 'None':
        # Apply correction and update p-values to the corrected ones
        decision, pvals, _, _ = multipletests(
            pvals,
            alpha=alpha,
            method=correction,
            is_sorted=False,
            returnsorted=False
        )
        # Return the rejected null hypothesis which have the highest
        diffs[np.where(decision == False)] = -np.Inf

    sorted_diff_indices = np.flip(np.argsort(diffs))
    last_positive = np.argmin(diffs[sorted_diff_indices] > 0)

    # accept only positive differences since those are markers
    # or until we reach max_n_genes

    # argmin updates only if it finds something smaller
    if last_positive == 0 and diffs[0] > 0:  # everything is positive
        last_positive = sorted_diff_indices.shape[0]
    final_indices = sorted_diff_indices[:min(last_positive, max_n_genes)]

    test_results = {'indices': final_indices,
                    'pvals': pvals[final_indices],
                    'diffs': diffs[final_indices]}
    return test_results
