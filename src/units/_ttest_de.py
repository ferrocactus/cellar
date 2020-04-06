import warnings

import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def _ttest_differential_expression(x_i, x_not_i, alpha=0.05, markers_n=200,
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

    markers_n: int
        Number of top markers to return. Note: return number could be smaller.

    correction: string
        Correction method to apply for multitest.
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    _______

    test_results: dict
        Dictionary of the form {
            'indices': array, shape (<=markers_n,)
            'pvals': array, shape (<=markers_n,)
            'diffs': array, shape (<=markers_n,)
        }
        where indices correspond to the indices of the columns which
        were found to be significant, pvals corresponds to the p-value
        as returned by multipletests(ttest), and diffs corresponds to
        the different in the means of the vectors x_in vs x_not_in
        for the given significant columns.

    """
    pvals = np.zeros(shape=(x_i.shape[1]))
    diffs = np.zeros_like(pvals)

    if x_i.shape[0] <= 1 or x_not_i.shape[0] <= 1:
        warnings.warn(
            "Array contains <=1 elements, ttest results may be inaccurate.")

    for j in range(x_i.shape[1]):  # current column index (gene)
        x_ij = x_i[:, j]
        x_not_ij = x_not_i[:, j]

        # WARNING: ttest will throw error if cluster has single point
        # Having equal_var=False runs Welch's t-test which assumes
        # different n and different variances
        _, pvals[j] = ttest_ind(x_ij, x_not_ij, equal_var=False)
        diffs[j] = np.mean(x_ij) - np.mean(x_not_ij)

    # Apply correction and update p-values to the corrected ones
    decision, pvals, _, _ = multipletests(
        pvals,
        alpha=alpha,
        method=correction,
        is_sorted=False,
        returnsorted=False
    )

    # Return the rejected null hypothesis which have the highest
    # difference of means of x_ij and x_not_ij
    diffs[diffs <= 0] = -np.Inf
    diffs[np.where(decision == False)] = -np.Inf
    sorted_diff_indices = np.flip(np.argsort(diffs))
    last_positive = np.argmin(diffs[sorted_diff_indices] > 0)

    # accept only positive differences since those are markers
    # or until we reach markers_n

    # argmin updates only if it finds something smaller
    if last_positive == 0 and diffs[0] > 0:  # everything is positive
        last_positive = sorted_diff_indices.shape[0]
    final_indices = sorted_diff_indices[:min(last_positive, markers_n)]

    test_results = {'indices': final_indices,
                    'pvals': pvals[final_indices],
                    'diffs': diffs[final_indices]}
    return test_results
