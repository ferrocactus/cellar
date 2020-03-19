import logging

import numpy as np
from sklearn.cluster import DBSCAN, KMeans, SpectralClustering

from ..methods import KMedoids
from ..utils.validation import _effective_n_clusters
from ._cluster_multiple import cluster_multiple
from ._unit import Unit

logger = logging.getLogger('Cluster')


def _get_wrapper(x, obj_def, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
    """
    Wrapper function for those classes which specify the number of clusters
    in advance and also have fit_predict implements. Classes include:
    KMedoids, KMeans, SpectralClustering.

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
    k, argtype = _effective_n_clusters(n_clusters)

    # If n_clusters determined to be single integer
    if argtype == 'int':
        y = obj_def(n_clusters=k, **kwargs).fit_predict(x)
        if eval_obj is not None:
            score = eval_obj.get(x, y)
            logger.info(f"Finished clustering with k={k}. Score={score:.2f}.")
        else:
            logger.info(f"Finished clustering with k={k}.")
        return y
    # If n_clusters determined to be a list of integers
    elif argtype == 'list':
        return cluster_multiple(
            x, obj_def=obj_def, k_list=k, attribute_name='n_clusters',
            eval_obj=eval_obj, method_name='fit_predict',
            n_jobs=n_jobs, **kwargs)


class Clu_KMedoids(Unit):
    """
    See src.methods._k_medoids
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
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
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Clusters and returns labels.

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
        logger.info("Initializing KMedoids.")
        return _get_wrapper(x, obj_def=KMedoids, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_KMeans(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
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
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Clusters and returns labels.

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
        logger.info("Initializing KMeans.")
        return _get_wrapper(x, obj_def=KMeans, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_SpectralClustering(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html
    """

    def __init__(self, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
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
        if 'affinity' not in kwargs:
            kwargs['affinity'] = 'nearest_neighbors'
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Clusters and returns labels.

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
        logger.info("Initializing SpectralClustering.")
        return _get_wrapper(x, obj_def=SpectralClustering, n_clusters=self.n_clusters,
                            eval_obj=self.eval_obj, n_jobs=self.n_jobs,
                            **self.kwargs)


class Clu_DBSCAN(Unit):
    """
    See https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """

    def __init__(self, eval_obj=None, n_jobs=None, **kwargs):
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
        self.eval_obj = eval_obj
        self.kwargs = kwargs

    def get(self, x):
        """
        Clusters and returns labels.

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
        logger.info("Initializing DBSCAN.")
        y = DBSCAN(**self.kwargs).fit_predict(x)
        unqy = len(np.unique(y))
        noise = np.sum(y == -1)
        if self.eval_obj is not None:
            score = self.eval_obj.get(x, y)
            logger.info(f"Found {unqy - (noise >= 1)} labels using DBSCAN."
                        "Score={score:.2f}.")
        else:
            logger.info(f"Found {unqy} labels using DBSCAN.")
        logger.info(f"Found {noise} noisy points. Assigning label -1.")
        return y
