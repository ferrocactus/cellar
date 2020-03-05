from ._unit import Unit

from abc import abstractmethod
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

ALPHA = 0.05
CORRECTION = 'hs'
MARKERS_N = 100


class Mark(Unit):
    """
    Base class for Marker Finding methods.
    """
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, **kwargs)
        self.name = 'Mark'

    @abstractmethod
    def get(self, x, labels):
        """
        Args:
            x (np.ndarray): 2D array
            labels (np.ndarray): 2D array of labels
        Returns:
        dict: {
            'label_1': {
                'indices': np.zeros(MARKERS_N),
                'pvals': np.zeros(MARKERS_N),
                'diffs': np.zeros(MARKERS_N)
            },
            'label_2': {
                ...
            }
            ...
        }
        """
        pass


class Mark_TTest(Mark):
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
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            alpha (float): Alpha value to use for hypothesis testing.
            correction (string): Correction method to apply for multitest.
            markers_n (int): Number of significant markers to return
                            (returned array could have less).
        """
        super().__init__(verbose, **kwargs)
        self.alpha = kwargs.get('alpha', ALPHA)
        self.correction = kwargs.get('correction', CORRECTION)
        self.markers_n = kwargs.get('markers_n', MARKERS_N)

    def get(self, x, labels, unq_labels=None):
        """
        Args:
            unq_labels (np.ndarray): np.unique(labels).
        """
        assert len(x.shape) == 2
        # Infer unique labels from input
        unq_labels = np.unique(labels) if unq_labels is None else unq_labels
        # To ensure non empty x_ij and x_not_ij
        if len(unq_labels) <= 1:
            self.vprint("Only one label found. Cannot run t-test.")
            return {}

        test_results = {}

        for i in unq_labels: # label to consider
            pvals = np.zeros(shape=(x.shape[1]))
            diffs = np.zeros_like(pvals)

            for j in range(x.shape[1]): # current column index (gene)
                x_ij = x[np.where(labels == i), j].squeeze() # Not empty
                x_not_ij = x[np.where(labels != i), j].squeeze() # Not empty

                #WARNING: ttest will throw error if cluster has single point
                # Having equal_var=False runs Welch's t-test which assumes
                # different n and different variances
                t, pvals[j] = ttest_ind(x_ij, x_not_ij, equal_var=False)
                diffs[j] = np.mean(x_ij) - np.mean(x_not_ij)

            # Apply correction and update p-values to the corrected ones
            decision, pvals, _, _ = multipletests(
                pvals,
                alpha=self.alpha,
                method=self.correction,
                is_sorted=False,
                returnsorted=False
            )

            # Return the rejected null hypothesis which have the highest
            # difference of means of x_ij and x_not_ij
            rejected_hs = np.where(decision == True)[0]
            sorted_diff_indices = np.flip(np.argsort(diffs[rejected_hs]))
            final_indices = rejected_hs[sorted_diff_indices][:self.markers_n]

            test_results[i] = {'indices': final_indices,
                               'pvals': pvals[final_indices],
                               'diffs': diffs[final_indices]}
            self.vprint("Finished finding markers for cluster {0}.".format(i))

        return test_results