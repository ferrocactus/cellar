import numpy as np
import Cluster_Ensembles as CE
from ..log import setup_logger
from ._unit import Unit
from ._cluster import Clu_KMeans
from ._cluster import Clu_KMedoids
from ._cluster import Clu_SpectralClustering
from ._cluster import Clu_Agglomerative
from ._cluster import Clu_DBSCAN
from ._cluster import Clu_Birch
from ._cluster import Clu_GaussianMixture
from ._cluster import Clu_Leiden
from ._cluster import Clu_Scanpy
from ..utils.validation import _validate_ensemble_methods


cluster_dict = {
    "KMeans": Clu_KMeans,
    "KMedoids": Clu_KMedoids,
    "Spectral": Clu_SpectralClustering,
    "Agglomerative": Clu_Agglomerative,
    "DBSCAN": Clu_DBSCAN,
    "Birch": Clu_Birch,
    "GaussianMixture": Clu_GaussianMixture,
    "Leiden": Clu_Leiden,
    "Scanpy": Clu_Scanpy
}


def clu_wrap(method):
    """
    Args:
        method (string): Method to use in the given step.
    Returns:
        object (Unit): Object of the right type.
    """
    if method not in cluster_dict:
        raise NotImplementedError("{0} method not implemented.".format(method))
    return cluster_dict[method]


class Ens_HyperGraph(Unit):
    """
    Ensemble clustering method based on hyper-graph partitioning.
    See https://github.com/GGiecold/Cluster_Ensembles
    """

    def __init__(self, methods=["KMedoids", "KMeans", "Spectral"],
                 n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        methods: list of clustering methods to use. Should be a list of strings.

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
        methods = _validate_ensemble_methods(methods)
        self.logger = setup_logger('Ensemble')
        if methods == "default":
            self.methods = ["KMedoids", "KMeans", "Spectral"]
        else:
            self.methods = methods
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
        self.logger.info("Initializing Ensemble Clustering.")
        self.logger.info("Using the following methods:")
        self.logger.info(", ".join(self.methods))

        if len(self.methods) == 1:
            # No need to do ensemble if only one method
            return clu_wrap(self.methods[0])(eval_obj=self.eval_obj,
                                    n_clusters=self.n_clusters,
                                    n_jobs=self.n_jobs, **self.kwargs).get(x)
        elif len(self.methods) < 1:
            raise ValueError("No methods specified for ensemble clustering.")

        # initialize empty partition matrix
        partitions = np.zeros((len(self.methods), x.shape[0]))

        for i, method in enumerate(self.methods):
            clu_obj = clu_wrap(method)(
                eval_obj=self.eval_obj, n_clusters=self.n_clusters,
                n_jobs=self.n_jobs, **self.kwargs
            )
            partitions[i, :] = clu_obj.get(x)

        ensemble_labels = CE.cluster_ensembles(partitions.astype(np.int))

        return ensemble_labels
