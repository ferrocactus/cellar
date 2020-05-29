import numpy as np
from joblib import Parallel, delayed

from ..log import setup_logger
from ..utils.validation import _validate_n_jobs
from ._ttest_de import _ttest_differential_expression
from ._unit import Unit


class DE_TTest(Unit):
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
        self.logger.info("Using {0} genes.".format(m))

        test_results = {}

        if _validate_n_jobs(self.n_jobs) == 1:
            for i in unq_labels:  # label to consider
                x_i = x[np.where(labels == i)]
                x_not_i = x[np.where(labels != i)]

                test_results[i] = _ttest_differential_expression(
                    x_i, x_not_i,
                    alpha=self.alpha,
                    markers_n=m,
                    correction=self.correction
                )
                self.logger.info(
                    "Finished finding DE genes for cluster={0}.".format(i))
        else:
            self.logger.info("Running DE gene discovery in parallel.")
            test_results = Parallel(n_jobs=self.n_jobs)(
                delayed(_ttest_differential_expression)(
                    x[np.where(labels == i)],
                    x[np.where(labels != i)],
                    alpha=self.alpha,
                    markers_n=m,
                    correction=self.correction
                ) for i in unq_labels
            )
            self.logger.info("Finished finding DE genes.")
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
