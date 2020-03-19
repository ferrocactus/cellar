import logging

import numpy as np
from sklearn.cluster import DBSCAN, KMeans, SpectralClustering

from ..methods import KMedoids
from ..utils.validation import _effective_n_clusters
from ._cluster_multiple import cluster_multiple
from ._unit import Unit

logger = logging.getLogger('Cluster')


class Clu_KMedoids(Unit):
    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        k, argtype = _effective_n_clusters(self.n_clusters)

        # If n_clusters determined to be single integer
        if argtype == 'int':
            y = KMedoids(n_clusters=k, **self.kwargs).fit_predict(x)
            score = eval_obj.get(x, y)
            logger.info(f"Finished clustering with k={k}. Score={score:.2f}.")
            return y
        # If n_clusters determined to be a list of integers
        elif argtype == 'list':
            return cluster_multiple(
                x, obj_def=KMedoids, k_list=k, attribute_name='n_clusters',
                eval_obj=self.eval_obj, method_name='fit_predict',
                n_jobs=self.n_jobs, **self.kwargs)


class Clu_KMeans(Unit):
    def __init__(self, **kwargs):
        super().__init__(verbose, name, **kwargs)
        self._obj = KMeans


class Clu_SpectralClustering(Unit):
    def __init__(self, **kwargs):
        super().__init__(verbose, name, **kwargs)
        if affinity not in self.kwargs:  # For consistency.
            self.kwargs['affinity'] = 'nearest_neighbors'
        self._obj = SpectralClustering


class Clu_DBSCAN(Unit):
    def __init__(self, **kwargs):
        super().__init__(verbose, name, **kwargs)
        self._obj = DBSCAN
