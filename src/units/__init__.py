from ._unit import Unit
from ._preprocess import Pre_Scanpy
from ._preprocess import Pre_ATAC
from ._dim_reduction import Dim_PCA
from ._dim_reduction import Dim_IncrementalPCA
from ._dim_reduction import Dim_KernelPCA
from ._dim_reduction import Dim_TruncatedSVD
from ._dim_reduction import Dim_UMAP
from ._dim_reduction import Dim_TSNE
from ._dim_reduction import Dim_MDS
from ._dim_reduction import Dim_FeatureAgglomeration
from ._dim_reduction import Dim_DiffusionMap
from ._dim_reduction import Dim_Isomap
from ._dim_reduction import Dim_SpectralEmbedding
from ._dim_reduction import Dim_UMAP_Paga
from ._cluster import Clu_KMeans
from ._cluster import Clu_KMedoids
from ._cluster import Clu_SpectralClustering
from ._cluster import Clu_Agglomerative
from ._cluster import Clu_DBSCAN
from ._cluster import Clu_Birch
from ._cluster import Clu_GaussianMixture
from ._cluster import Clu_Leiden
from ._cluster_ensemble import Ens_HyperGraph
from ._evaluation import Eval_Silhouette
from ._evaluation import Eval_DaviesBouldin
from ._evaluation import Eval_CalinskiHarabasz
from ._de import DE_TTest_Cellar
from ._de import DE_TTest
from ._de import DE_Rank
from ._de import DE_LRT
from ._de import DE_Wald
from ._converter import Con
from ._converter import convert
from ._identificator import Ide_HyperGeom
from ._ss_cluster import SSClu_SeededKMeans
from ._ss_cluster import SSClu_ConstrainedKMeans
from ._align import Ali_Scanpy_Ingest

from ..utils.exceptions import MethodNotImplementedError

translation_dict = {
    "preprocess": {
        "Scanpy": Pre_Scanpy,
        "sc-ATAC-seq": Pre_ATAC
    },
    "dim_reduction": {
        "PCA": Dim_PCA,
        "Incremental PCA": Dim_IncrementalPCA,
        "Kernel PCA": Dim_KernelPCA,
        "Truncated SVD": Dim_TruncatedSVD,
        "UMAP": Dim_UMAP,
        "UMAP + Paga": Dim_UMAP_Paga,
        "MDS": Dim_MDS,
        "Feature Agglomeration": Dim_FeatureAgglomeration,
        "Diffusion Map": Dim_DiffusionMap,
        "Isomap": Dim_Isomap,
        "Spectral Embedding": Dim_SpectralEmbedding
        # "Autoencoder": Dim_AE
    },
    "cluster": {
        "KMeans": Clu_KMeans,
        "KMedoids": Clu_KMedoids,
        "Spectral": Clu_SpectralClustering,
        "Agglomerative": Clu_Agglomerative,
        # "DBSCAN": Clu_DBSCAN,
        # "Birch": Clu_Birch,
        "GaussianMixture": Clu_GaussianMixture,
        "Leiden": Clu_Leiden,
        "Ensemble": Ens_HyperGraph
    },
    "cluster_eval": {
        "Silhouette": Eval_Silhouette,
        "DaviesBouldin": Eval_DaviesBouldin,
        "CalinskiHarabasz": Eval_CalinskiHarabasz
    },
    "de": {
        "Cellar (TTest)": DE_TTest_Cellar,
        "TTest": DE_TTest,
        "Rank": DE_Rank,
        "Likelihood Ratio": DE_LRT,
        "Wald": DE_Wald
    },
    "conversion": {
        "Converter": Con
    },
    "identification": {
        "HyperGeom": Ide_HyperGeom
    },
    "ss_cluster": {
        "SeededKMeans": SSClu_SeededKMeans,
        "ConstrainedKMeans": SSClu_ConstrainedKMeans
    },
    "visualization": {
        "UMAP": Dim_UMAP,
        "UMAP + Paga": Dim_UMAP_Paga,
        "TSNE": Dim_TSNE,
        "PCA": Dim_PCA,
        "Kernel PCA": Dim_KernelPCA,
        "Feature Agglomeration": Dim_FeatureAgglomeration,
        "MDS": Dim_MDS,
        "Diffusion Map": Dim_DiffusionMap,
        "Isomap": Dim_Isomap,
        "Spectral Embedding": Dim_SpectralEmbedding
    },
    "align": {
        "Scanpy Ingest": Ali_Scanpy_Ingest
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
    'Pre_Scanpy',
    'Pre_ATAC',
    'Dim_PCA',
    'Dim_IncrementalPCA',
    'Dim_KernelPCA',
    'Dim_TruncatedSVD',
    'Dim_UMAP',
    'Dim_UMAP_Paga',
    'Dim_TSNE',
    'Dim_MDS',
    'Dim_Isomap',
    'Dim_FeatureAgglomeration',
    'Dim_DiffusionMap',
    'Dim_SpectralEmbedding',
    'Clu_KMeans',
    'Clu_KMedoids',
    'Clu_SpectralClustering',
    'Clu_Agglomerative',
    # 'Clu_DBSCAN',
    # 'Clu_Birch',
    'Clu_GaussianMixture',
    'Clu_Leiden',
    'Ens_HyperGraph',
    'Eval_Silhouette',
    'Eval_DaviesBouldin',
    'Eval_CalinskiHarabasz',
    'DE_TTest',
    'DE_Rank',
    'DE_LRT',
    'DE_Wald',
    'Con',
    'Ide_HyperGeom',
    'SSClu_SeededKMeans',
    'SSClu_ConstrainedKMeans',
    'Ali_Scanpy_Ingest',
    'convert'
]
