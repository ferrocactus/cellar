from ._unit import Unit
from ._cluster import (Clu_KMeans, Clu_KMedoids, Clu_SpectralClustering,
                       Clu_Agglomerative, Clu_DBSCAN, Clu_Birch,
                       Clu_GaussianMixture)
from ._dim_reduction import Dim_PCA, Dim_UMAP, Dim_TSNE
from ._evaluation import (Eval_Silhouette, Eval_DaviesBouldin,
                          Eval_CalinskiHarabasz)
from ._markers import Mark_TTest
from ._converter import Con
from ._identificator import Ide_HyperGeom
from ._ss_cluster import SSClu_SeededKMeans
from ._ss_cluster import SSClu_ConstrainedKMeans
from ._ss_cluster import SSClu_ConstrainedSeededKMeans

translation_dict = {
    "cluster": {
        "KMeans": Clu_KMeans,
        "KMedoids": Clu_KMedoids,
        "Spectral": Clu_SpectralClustering,
        "Agglomerative": Clu_Agglomerative,
        "DBSCAN": Clu_DBSCAN,
        "Birch": Clu_Birch,
        "GaussianMixture": Clu_GaussianMixture
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
    "markers": {
        "TTest": Mark_TTest
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
    }
}

__all__ = [
    'translation_dict',
    'Clu_KMeans',
    'Clu_KMedoids',
    'Clu_SpectralClustering',
    'Clu_Agglomerative',
    'Clu_DBSCAN',
    'Clu_Birch',
    'Clu_GaussianMixture',
    'Dim_PCA',
    'Dim_UMAP',
    'Dim_TSNE',
    'Eval_Silhouette',
    'Eval_DaviesBouldin',
    'Eval_CalinskiHarabasz',
    'Mark_TTest',
    'Con',
    'Ide_HyperGeom',
    'SSClu_SeededKMeans',
    'SSClu_ConstrainedKMeans',
    'SSClu_ConstrainedSeededKMeans'
]
