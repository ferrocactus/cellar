from abc import abstractmethod

import numpy as np
from sklearn.cluster import KMeans
from umap import UMAP

from ..methods._kmeans import KMeans as ConstrainedKMeans
from ..log import setup_logger
from ._unit import Unit

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

    def get(self, x, labels, mask):
        """
        Find the non-negative values and let them define the n_clusters.
        Args:
            x (np.ndarray): Data points (n_points x n_features)
            labels (np.ndarray): Labels (n_points x 1)
            mask: ignored, used for consistencyd
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
    def __init__(self, **kwargs):
        self.logger = setup_logger('ConstrainedKMeans')
        self.kwargs = kwargs

    def get(self, x, labels, mask):
        unq_labels = np.unique(labels)
        n_clusters = len(unq_labels)

        mask = mask.astype(np.int32)
        labels = labels.astype(np.int32)

        constrainedkmeans = ConstrainedKMeans(n_clusters=n_clusters,
                                    mask=mask, fixed_labels=labels)
        labels = constrainedkmeans.fit_predict(x)
        return labels


class SSClu_ConstrainedSeededKMeans(Unit):
    def __init__(self, **kwargs):
        self.logger = setup_logger('ConstrainedSeededKMeans')
        self.kwargs = kwargs

    def get(self, x, labels, mask):
        unq_labels = np.unique(labels)
        n_clusters = len(unq_labels)
        centroids = []
        self.logger.info(
            "Found {0} unique labels. Using constrained seeded KMeans.".format(n_clusters))

        for i in range(n_clusters):
            centroid = np.mean(x[labels == unq_labels[i]], axis=0)
            centroids.append(centroid)

        centroids = np.array(centroids)
        mask = mask.astype(np.int32)
        labels = labels.astype(np.int32)

        constrainedkmeans = ConstrainedKMeans(n_clusters=n_clusters,
                                    init=centroids, n_init=1,
                                    mask=mask, fixed_labels=labels)
        labels = constrainedkmeans.fit_predict(x)
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
