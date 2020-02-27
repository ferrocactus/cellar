from ._cluster import (
    Clu_KMedoids,
    Clu_KMeans,
    Clu_SpectralClustering
)
from ._dim_reduction import (
    Dim_PCA,
    Dim_UMAP,
    Dim_TSNE
)
from ._evaluation import (
    Eval_SilhouetteScore,
    Eval_DaviesBouldinScore
)
from ._markers import (
    Mark_TTest
)
from ._converter import (
    Con
)

translation_dict = {
    "cluster": {
        "KMedoids": Clu_KMedoids,
        "KMeans": Clu_KMeans,
        "SpectralClustering": Clu_SpectralClustering
    },
    "dim_reduction": {
        "PCA": Dim_PCA,
        "UMAP": Dim_UMAP,
        "TSNE": Dim_TSNE
    },
    "cluster_eval": {
        "SilhouetteScore": Eval_SilhouetteScore,
        "DaviesBouldinScore": Eval_DaviesBouldinScore
    },
    "markers": {
        "TTest": Mark_TTest
    },
    "converter": {
        "Converter": Con
    }
}