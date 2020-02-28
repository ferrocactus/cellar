from ._unit import Unit

from abc import abstractmethod

#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans
#from copkmeans.cop_kmeans import cop_kmeans

class SSClu(Unit):
    """
    Base class for Semi-Supervised Clustering.
    """
    def __init__(self, verbose=False, **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, **kwargs)
        self.name = 'Dim'

    @abstractmethod
    def get(self, x):
        """
        Returns the labels of x.

        Args:
            x (np.ndarray): Data in matrix (n x d) form.
        Returns:
            (np.ndarray): The labels of x.
        """
        pass


class SSClu_COPKMeans(SSClu):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)

    def get(self, x, n, ml, cl):
        clusters, _ = cop_kmeans(x, k=n, ml=ml, cl=cl)
        return clusters


class SSClu_PCKMeans(SSClu):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)

    def get(self, x, n, ml, cl):
        pckmeans = PCKMeans(n_clusters=n)
        pckmeans.fit(x, ml=ml, cl=cl)
        return pckmeans.labels_