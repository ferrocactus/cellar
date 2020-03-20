import warnings
from abc import abstractmethod

import numpy as np
from joblib import Parallel, delayed
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from ..log import setup_logger
from ._unit import Unit


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
        warnings.warn("Array contains <=1 elements, ttest results may be inaccurate.")

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
    if last_positive == 0 and diffs[0] > 0:
        last_positive = sorted_diff_indices.shape[0]
    final_indices = sorted_diff_indices[:min(last_positive, markers_n)]

    test_results = {'indices': final_indices,
                    'pvals': pvals[final_indices],
                    'diffs': diffs[final_indices]}
    return test_results


class Mark_TTest(Unit):
    """
    One-vs-all T-Test.

    Two sided test for the null H_0: The two populations have the same average.
    If the p-value is small (<alpha), reject the null. Apply correction for
    multiple testing.

    Takes a two dimensional matrix x where columns correspond to genes and
    labels for every cell, i.e., row. For every cluster label, iterate through
    all columns (genes) and isolate the rows (cells) which are assigned that
    particular label. Then run Welch's T-Test to determine if this particular
    gene is any significant for the current label. In addition to requiring
    the p-value being less than alpha (corrected for multiple testing) we
    choose the top_k genes which exhibit the highest difference of the means.
    """

    def __init__(self, alpha=0.05, markers_n=200, correction='hs', n_jobs=None):
        """
        Parameters
        __________

        alpha: float, between 0 and 1
        Alpha value to use for hypothesis testing.

        markers_n: int
            Number of top markers to return. Note: return number could be smaller.

        correction: string
            Correction method to apply for multitest.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.
            If None or 1, will not run multithreading.

        """
        self.logger = setup_logger('TTest')
        self.alpha = alpha
        self.markers_n = markers_n
        self.correction = correction
        self.n_jobs = n_jobs

    def get(self, x, labels, unq_labels=None):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        labels: array, shape (n_samples,)
            Labels of data points.

        unq_labels: array or None, shape (np.unique(labels),)
            Will get recomputed, but can pass as argument for performance.

        Returns
        _______

        test_results: dict
            Dictionary of the form {
                label {
                    'indices': array, shape (<=markers_n,)
                    'pvals': array, shape (<=markers_n,)
                    'diffs': array, shape (<=markers_n,)
                },
                ...
            }
        """
        # Infer unique labels from input
        unq_labels = np.unique(labels) if unq_labels is None else unq_labels
        # To ensure non empty x_ij and x_not_ij
        if len(unq_labels) <= 1:
            self.logger.info("Only one label found. Cannot run t-test.")
            return {}

        m = int(self.markers_n*x.shape[1]
                ) if self.markers_n < 1 else self.markers_n
        self.logger.info(f"Using {m} markers.")

        test_results = {}

        if self.n_jobs is None or self.n_jobs == 1:
            for i in unq_labels:  # label to consider
                x_i = x[np.where(labels == i)]
                x_not_i = x[np.where(labels != i)]

                test_results[i] = _ttest_differential_expression(
                    x_i, x_not_i,
                    alpha=self.alpha,
                    markers_n=m,
                    correction=self.correction
                )
                self.logger.info(f"Finished finding markers for cluster={i}.")
        else:
            self.logger.info("Running marker discovery in parallel.")
            test_results = Parallel(n_jobs=self.n_jobs)(
                delayed(_ttest_differential_expression)(
                    x[np.where(labels == i)],
                    x[np.where(labels != i)],
                    alpha=self.alpha,
                    markers_n=m,
                    correction=self.correction
                ) for i in unq_labels
            )
            self.logger.info("Finished finding markers.")
            test_results = dict(zip(unq_labels, test_results))

        return test_results

    def get_subset(self, x, subset_indices):
        """
        Do differential expression for the points specified in subset_indices
        versus the rest of the population.
        """
        m = int(self.markers_n*x.shape[1]
                ) if self.markers_n < 1 else self.markers_n
        mask = np.zeros(len(x), dtype=bool)
        mask[subset_indices, ] = True
        return {'0': _ttest_differential_expression(
            x[mask],
            x[~mask],
            alpha=self.alpha,
            markers_n=m,
            correction=self.correction
        )}
