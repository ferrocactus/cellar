from ._cluster import *
from ._dim_reduction import *
from ._evaluation import *
from ._markers import *
from ._converter import *
from ._identificator import *
from ._ss_cluster import *

translation_dict = {
    "cluster": {
        "KMedoids": Clu_KMedoids,
        "KMeans": Clu_KMeans,
        "SpectralClustering": Clu_SpectralClustering
    },
    "dim_reduction": {
        "PCA": Dim_PCA,
        "UMAP": Dim_UMAP,
        "TSNE": Dim_TSNE,
        "Autoencoder": Dim_AE
    },
    "cluster_eval": {
        "SilhouetteScore": Eval_SilhouetteScore,
        "DaviesBouldinScore": Eval_DaviesBouldinScore
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
        "UMAP": SSClu_UMAP,
        "COPKMeans": SSClu_COPKMeans,
        "PCKMeans": SSClu_PCKMeans
    }
}