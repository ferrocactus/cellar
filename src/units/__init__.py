from ._unit import Unit
from ._cluster import (Clu_KMeans, Clu_KMedoids, Clu_SpectralClustering,
                       Clu_Agglomerative, Clu_DBSCAN, Clu_Birch)
from ._dim_reduction import Dim_PCA, Dim_UMAP, Dim_TSNE
from ._evaluation import (Eval_SilhouetteScore, Eval_DaviesBouldinScore,
                          Eval_CalinskiHarabasz)
from ._markers import Mark_TTest
from ._converter import Con
from ._identificator import Ide_HyperGeom
from ._ss_cluster import SSClu_SeededKMeans

translation_dict = {
    "cluster": {
        "KMeans": Clu_KMeans,
        "KMedoids": Clu_KMedoids,
        "SpectralClustering": Clu_SpectralClustering,
        "Agglomerative": Clu_Agglomerative,
        "DBSCAN": Clu_DBSCAN,
        "Birch": Clu_Birch
    },
    "dim_reduction": {
        "PCA": Dim_PCA,
        "UMAP": Dim_UMAP,
        "TSNE": Dim_TSNE
        # "Autoencoder": Dim_AE
    },
    "cluster_eval": {
        "SilhouetteScore": Eval_SilhouetteScore,
        "DaviesBouldinScore": Eval_DaviesBouldinScore,
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
        "SeededKMeans": SSClu_SeededKMeans
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
    'Dim_PCA',
    'Dim_UMAP',
    'Dim_TSNE',
    'Eval_SilhouetteScore',
    'Eval_DaviesBouldinScore',
    'Eval_CalinskiHarabasz',
    'Mark_TTest',
    'Con',
    'Ide_HyperGeom',
    'SSClu_SeededKMeans'
]
