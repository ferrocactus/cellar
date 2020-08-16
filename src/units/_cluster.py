import logging

import anndata
from anndata import AnnData
import leidenalg
import igraph
import numpy as np
import scanpy
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.cluster import Birch
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import kneighbors_graph

from ..log import setup_logger
from ..methods import KMedoids
from ..utils.validation import _validate_clu_n_clusters
from ._cluster_multiple import cluster_multiple
from ._unit import Unit
from ._evaluation import Eval_Silhouette


default_eval_obj = Eval_Silhouette()


def _get_wrapper(x, obj_def, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj,
                 n_jobs=None, attribute_name='n_clusters', **kwargs):
    """
    Wrapper function for those classes which specify the number of clusters
    in advance and also have fit_predict implemented. Classes include:
    KMedoids, KMeans, SpectralClustering, AgglomerativeClustering, Birch.

    Parameters
    __________

    x: array, shape (n_samples, n_features)
        The data array.

    obj_def: object name
        Object to be instantiated in this function.

    n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
        Array containing the different values of clusters to try,
        or single int specifying the number of clusters,
        or tuple of the form (a, b, c) which specifies a range
        for (x=a; x<b; x+=c)

    eval_obj: Eval or None, default None
        Evaluation object to compare performance of different trials.

    n_jobs: int or None, default None
        Number of jobs to use if multithreading. See
        https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

    attribute_name: string, default 'n_clusters'
        Name of the obj.attribute_name that corresponds to n_clusters.

    **kwargs: dictionary
        Dictionary of parameters that will get passed to obj_def
        when instantiating it.

    Returns
    _______

    y: array, shape (n_samples,)
        List of labels that correspond to the best clustering k, as
        evaluated by eval_obj.

    """
    # Determine type of n_clusters passed
    k = _validate_clu_n_clusters(n_clusters, x.shape[0])

    # If n_clusters determined to be single integer
    if isinstance(k, int):
        logger = setup_logger('Cluster.Single')
        kwargs[attribute_name] = k
        y = obj_def(**kwargs).fit_predict(x)
        if eval_obj is None:
            eval_obj = default_eval_obj
        score = eval_obj.get(x, y)
        logger.info(
            "Finished clustering with k={0}. Score={1:.2f}.".format(k,
                                                                    score))
        return y, score
    # If n_clusters determined to be a list of integers
    elif isinstance(k, (list, np.ndarray)):
        return cluster_multiple(
            x, obj_def=obj_def, k_list=k, attribute_name=attribute_name,
            eval_obj=eval_obj, method_name='fit_predict',
            n_jobs=n_jobs, **kwargs)


class Clu_KMedoids(Unit):
    """
    See src.methods._k_medoids
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('KMedoids')
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing KMedoids.")
        return _get_wrapper(x, obj_def=KMedoids, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_KMeans(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('KMeans')
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing KMeans.")
        return _get_wrapper(x, obj_def=KMeans, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_SpectralClustering(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('Spectral Clustering')
        if 'affinity' not in kwargs:
            kwargs['affinity'] = 'nearest_neighbors'
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing SpectralClustering.")
        return _get_wrapper(x, obj_def=SpectralClustering,
                            n_clusters=self.n_clusters, eval_obj=self.eval_obj,
                            n_jobs=self.n_jobs, **self.kwargs)


class Clu_Agglomerative(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('Agglomerative')
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing Agglomerative Clustering.")
        return _get_wrapper(x, obj_def=AgglomerativeClustering,
                            n_clusters=self.n_clusters, eval_obj=self.eval_obj,
                            n_jobs=self.n_jobs, **self.kwargs)


class Clu_Birch(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.Birch.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('Birch')
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing Birch Clustering.")
        return _get_wrapper(x, obj_def=Birch, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_DBSCAN(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """

    def __init__(self, eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        eval_obj: Eval or None, default None
            Evaluation object to evaluate clustering.

        n_jobs: Ignored

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('DBSCAN')
        self.eval_obj = eval_obj
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing DBSCAN.")
        y = DBSCAN(**self.kwargs).fit_predict(x)
        unqy = len(np.unique(y))
        noise = np.sum(y == -1)
        if self.eval_obj is not None:
            score = self.eval_obj.get(x, y)
            self.logger.info(
                "Found {0} labels using DBSCAN.".format(unqy-(noise >= 1)) +
                "Score={0:.2f}.".format(score))
        else:
            self.logger.info("Found {0} labels using DBSCAN.".format(unqy))
        self.logger.info(
            "Found {0} noisy points. Assigning label -1.".format(noise))
        return y, score


class Clu_GaussianMixture(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=default_eval_obj, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('GaussianMixture')
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing Gaussian Mixture Model.")
        return _get_wrapper(x, obj_def=GaussianMixture, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            attribute_name='n_components', **self.kwargs)


class Clu_Leiden(Unit):
    """
    See https://github.com/vtraag/leidenalg
    """

    def __init__(self, n_neighbors=10, eval_obj=default_eval_obj, **kwargs):
        """
        Parameters
        __________
        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        self.logger = setup_logger('Leiden')
        self.n_neighbors = n_neighbors
        self.eval_obj = eval_obj
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        self.logger.info("Initializing Leiden Clustering.")
        try:
            import warnings
            from numba.errors import NumbaPerformanceWarning
            warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
        except:
            pass

        is_anndata = isinstance(x, AnnData)

        if not is_anndata:
            adata = anndata.AnnData(X=x)
            # Use 0 pcs since x is supposed to be an embedding already
            n_pcs = 0
        else:
            adata = x # reference
            adata.obsm['X_pca'] = adata.obsm['x_emb'].copy()
            n_pcs = None

        scanpy.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
        scanpy.tl.leiden(adata)
        labels = np.squeeze(np.array(adata.obs['leiden'])).astype(np.int)

        if is_anndata:
            adata.obsm.pop('X_pca')
        return labels, 0


# class Clu_Scanpy(Unit):
#     """
#     See https://scanpy.readthedocs.io
#     """

#     def __init__(self, **kwargs):
#         """
#         Parameters
#         __________
#         **kwargs: dictionary
#             Dictionary of parameters that will get passed to obj_def
#             when instantiating it.

#         """
#         self.logger = setup_logger('Scanpy')
#         self.kwargs = kwargs

#     def get(self, x):
#         """
#         Parameters
#         __________
#         x: array, shape (n_samples, n_features)
#             The data array.

#         Returns
#         _______
#         y: array, shape (n_samples,)
#             List of labels that correspond to the best clustering k, as
#             evaluated by eval_obj.

#         """

#         self.logger.info("Initializing Scanpy.")
#         ann = anndata.AnnData(X=x)
#         scanpy.pp.neighbors(ann, n_neighbors=10, n_pcs=40)
#         scanpy.tl.leiden(ann)
#         labels = np.squeeze(np.array(ann.obs)).astype(np.int)
#         return labels, self.eval_obj.get(x, labels)
