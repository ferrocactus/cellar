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
from ._cluster_ensemble import Ens_HyperGraph
from ._dim_reduction import Dim_PCA
from ._dim_reduction import Dim_UMAP
from ._dim_reduction import Dim_TSNE
from ._evaluation import Eval_Silhouette
from ._evaluation import Eval_DaviesBouldin
from ._evaluation import Eval_CalinskiHarabasz
from ._de import DE_TTest
from ._converter import Con
from ._identificator import Ide_HyperGeom
from ._ss_cluster import SSClu_SeededKMeans
from ._ss_cluster import SSClu_ConstrainedKMeans
from ._ss_cluster import SSClu_ConstrainedSeededKMeans

from ..utils.exceptions import MethodNotImplementedError

translation_dict = {
    "cluster": {
        "KMeans": Clu_KMeans,
        "KMedoids": Clu_KMedoids,
        "Spectral": Clu_SpectralClustering,
        "Agglomerative": Clu_Agglomerative,
        # "DBSCAN": Clu_DBSCAN,
        # "Birch": Clu_Birch,
        "GaussianMixture": Clu_GaussianMixture,
        "Leiden": Clu_Leiden,
        "Scanpy": Clu_Scanpy,
        "Ensemble": Ens_HyperGraph
    },
    "dim_reduction": {
        "PCA": Dim_PCA,
        "UMAP": Dim_UMAP,
        "TSNE": Dim_TSNE
        # "Autoencoder": Dim_AE
    },
    "cluster_eval": {
        "Silhouette": Eval_Silhouette,
        "DaviesBouldin": Eval_DaviesBouldin,
        "CalinskiHarabasz": Eval_CalinskiHarabasz
    },
    "de": {
        "TTest": DE_TTest
    },
    "conversion": {
        "Converter": Con
    },
    "identification": {
        "HyperGeom": Ide_HyperGeom
    },
    "ss_cluster": {
        "SeededKMeans": SSClu_SeededKMeans,
        "ConstrainedKMeans": SSClu_ConstrainedKMeans,
        "ConstrainedSeededKMeans": SSClu_ConstrainedSeededKMeans
    },
    "visualization": {
        "UMAP": Dim_UMAP,
        "TSNE": Dim_TSNE
    }
}


def _method_exists(step, method):
    """
    Check if method in step has been implemented.
    """
    if step not in translation_dict:
        raise MethodNotImplementedError(
            "{0} step not implemented.".format(step))
    if method not in translation_dict[step]:
        raise MethodNotImplementedError(
            "{0} method not implemented.".format(method))


def wrap(step, method):
    """
    Wrapper function that takes a pipeline step and method
    and returns the corresponding object.

    Args:
        step (string): The step in the pipeline.
        method (string): Method to use in the given step.
    Returns:
        object (Unit): Object of the right type.
    """
    try:
        _method_exists(step, method)
    except MethodNotImplementedError as e:
        return str(e)
    return translation_dict[step][method]


__all__ = [
    'translation_dict',
    '_method_exists',
    'wrap',
    'Clu_KMeans',
    'Clu_KMedoids',
    'Clu_SpectralClustering',
    'Clu_Agglomerative',
    # 'Clu_DBSCAN',
    # 'Clu_Birch',
    'Clu_GaussianMixture',
    'Clu_Leiden',
    'Clu_Scanpy',
    'Ens_HyperGraph',
    'Dim_PCA',
    'Dim_UMAP',
    'Dim_TSNE',
    'Eval_Silhouette',
    'Eval_DaviesBouldin',
    'Eval_CalinskiHarabasz',
    'DE_TTest',
    'Con',
    'Ide_HyperGeom',
    'SSClu_SeededKMeans',
    'SSClu_ConstrainedKMeans',
    'SSClu_ConstrainedSeededKMeans'
]
