from ._unit import Unit

from abc import abstractmethod

from umap import UMAP
import numpy as np
#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans
#from copkmeans.cop_kmeans import cop_kmeans

class SSClu(Unit):
    """
    Base class for Semi-Supervised Clustering.
    """
    def __init__(self, verbose=False, name='SSClu', **kwargs):
        """
        Args:
            verbose (bool): Printing flag.
            **kwargs: Argument dict.
        """
        super().__init__(verbose, name, **kwargs)

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


class SSClu_UMAP(SSClu):
    def __init__(self, verbose=False, name='SS UMAP', **kwargs):
        super().__init__(verbose, name, **kwargs)

    def get(self, x, labels, clu, eval):
        umap = UMAP(**self.kwargs)
        self.vprint("Finding embeddings.")
        emb = umap.fit_transform(x, y=labels)
        new_labels = clu.get(emb, eval)
        ind = np.where(labels != -1)
        new_labels[ind] = labels[ind]
        return emb, new_labels


class SSClu_COPKMeans(SSClu):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)

    def get(self, x, n, ml, cl):
        pass
        clusters, _ = cop_kmeans(x, k=n, ml=ml, cl=cl)
        return clusters


class SSClu_PCKMeans(SSClu):
    def __init__(self, verbose=False, **kwargs):
        super().__init__(verbose, **kwargs)

    def get(self, x, n, ml, cl):
        pass
        pckmeans = PCKMeans(n_clusters=n)
        pckmeans.fit(x, ml=ml, cl=cl)
        return pckmeans.labels_