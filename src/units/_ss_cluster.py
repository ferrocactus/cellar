from abc import abstractmethod

import numpy as np
from sklearn.cluster import KMeans
from umap import UMAP

from ..methods import ConstrainedKMeans
from ..log import setup_logger
from ._unit import Unit
from ..utils.exceptions import InappropriateArgument
from ..utils.validation import _validate_n_clusters
from ..utils.validation import _validate_cluster_list


class SSClu_SeededKMeans(Unit):
    """
    Similar to KMeans, given a dataset x, generate e K-partitioning
    of x so that the KMeans objective is locally minimized. The difference
    lies in the way we initialize centroids. Given a list of labels for every
    point in x, we let the number of clusters equal the number of unique
    non-negative labels, and for every such cluster, we initialize a centroid
    based on the average of the points with that label.

    Source: Basu, Sugato, et al. “Semi-Supervised Clustering by Seeding.”
    Proceedings of the 19th International Conference on Machine Learning
    (ICML-2002), no. July, 2002, pp. 19–26.
    """

    def __init__(self, **kwargs):
        self.logger = setup_logger('SeededKMeans')
        self.kwargs = kwargs

    def get(self, x, labels, preserved_labels=None):
        """
        Find the non-negative values and let them define the n_clusters.
        Args:
            x (np.ndarray): Data points (n_points x n_features)
            labels (np.ndarray): Labels (n_points x 1)
            preserved_labels: ignored, used for consistency
        Returns:
            (np.ndarray): New labels after clustering (n_points x 1)
        """
        unq_labels = np.unique(labels[labels >= 0])
        n_clusters = len(unq_labels)
        centroids = []
        self.logger.info(
            "Found {0} unique labels. Using seeded KMeans.".format(n_clusters))

        for i in range(n_clusters):
            centroid = np.mean(x[labels == unq_labels[i]], axis=0)
            centroids.append(centroid)

        centroids = np.array(centroids)

        kmeans = KMeans(n_clusters=n_clusters, init=centroids, n_init=1)
        labels = kmeans.fit_predict(x)
        return labels


class SSClu_ConstrainedKMeans(Unit):
    """
    See documentation in src/methods/ConstrainedKMeans.py.
    """

    def __init__(self, n_clusters, **kwargs):
        self.logger = setup_logger('Constrained KMeans')
        self.n_clusters = int(n_clusters)
        self.kwargs = kwargs

    def get(self, x, labels, preserved_labels):
        unq_labels = np.unique(labels)
        _validate_n_clusters(self.n_clusters, x.shape[0])
        if (self.n_clusters < len(preserved_labels)):
            raise InappropriateArgument("Number of clusters is less ",
                                        "than the number of preserved labels")

        can_change = (1 - np.isin(labels, preserved_labels)).astype(np.int32)
        init_labels = np.asarray(labels).astype(np.int32)

        ckm = ConstrainedKMeans(n_clusters=self.n_clusters, **self.kwargs)
        labels = ckm.fit_predict(
            x, can_change=can_change, init_labels=init_labels)
        return labels


class SSClu_UMAP(Unit):
    def __init__(self, **kwargs):
        self.logger = setup_logger('UMAP')
        self.kwargs = kwargs

    def get(self, x, labels, clu, eval):
        umap = UMAP(**self.kwargs)
        self.logger.info("Finding embeddings.")
        emb = umap.fit_transform(x, y=labels)
        new_labels = clu.get(emb, eval)
        ind = np.where(labels != -1)
        new_labels[ind] = labels[ind]
        return new_labels
