from abc import abstractmethod

import numpy as np
from sklearn.cluster import KMeans
from umap import UMAP

from ._unit import Unit

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


class SSClu_SeededKMeans(SSClu):
    """
    Similar to KMeans, given a dataset x, generate e K-partitioning
    of x so that the KMeans objective is locally minimized. The difference
    lies in the way we initialize centroids. Given a list of labels for every
    point in x, we let the number of clusters equal the number of unique
    non-negative labels, and for every such cluster, we initialize a centroid
    based on the average of the points with that label.

    Source: Basu, Sugato, et al. “Semi-Supervised Clustering by Seeding.”
    Proceedings of the 19th International Conference on Machine Learning (ICML-2002),
    no. July, 2002, pp. 19–26.
    """

    def __init__(self, **kwargs):
        self.name = 'SeededKMeans'
        self.kwargs = kwargs

    def get(self, x, labels):
        """
        Find the non-negative values and let them define the n_clusters.
        Args:
            x (np.ndarray): Data points (n_points x n_features)
            labels (np.ndarray): Labels (n_points x 1)
        Returns:
            (np.ndarray): New labels after clustering (n_points x 1)
        """
        unq_labels = np.unique(labels[labels >= 0])
        n_clusters = len(unq_labels)
        centroids = []
        self.logger.info(f"Found {n_clusters} unique labels. Using seeded KMeans.")

        for i in range(n_clusters):
            centroid = np.mean(x[labels == unq_labels[i]], axis=0)
            centroids.append(centroid)

        centroids = np.array(centroids)

        kmeans = KMeans(n_clusters=n_clusters, init=centroids, n_init=1)
        labels = kmeans.fit_predict(x)
        return labels


class SSClu_UMAP(SSClu):
    def __init__(self, **kwargs):
        self.name = 'SS_UMAP'
        self.kwargs = kwargs

    def get(self, x, labels, clu, eval):
        umap = UMAP(**self.kwargs)
        self.logger.info("Finding embeddings.")
        emb = umap.fit_transform(x, y=labels)
        new_labels = clu.get(emb, eval)
        ind = np.where(labels != -1)
        new_labels[ind] = labels[ind]
        return new_labels


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
