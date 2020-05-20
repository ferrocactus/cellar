defaults <- list(
    "dim" = '10',
    "clu" = '(3, 5, 1)',
    "mark_alpha" = 0.05,
    "mark_no" = 200
)

options <- list(
    "dim" = c("PCA", "UMAP", "TSNE"),
    "clu" = c("KMeans", "KMedoids", "Spectral", "Agglomerative",
              "GaussianMixture", "Leiden", "Scanpy", "Ensemble"),
    "clu_ensemble" = c("All"="All", "KMeans"="KMeans", "KMedoids"="KMedoids",
                       "GaussianMixture"="GaussianMixture",
                       "Spectral"="Spectral", "Agglomerative"="Agglomerative",
                       "Leiden"="Leiden", "Scanpy"="Scanpy"),
    "clu_no_n_clusters" = c("DBSCAN", "Leiden", "Scanpy"),
    "eval" = c("Silhouette", "DaviesBouldin", "CalinskiHarabasz"),
    "correction" = c("holm-sidak", "bonferroni", "sidak", "holm",
                     "simes-hochberg", "hommel", "fdr_bh", "fdr_by",
                     "fdr_tsbh", "fdr_tsbky"),
    "converter" = c("id-to-name", "name-to-id"),
    "tissues" = c("all", "Spleen", "Thyroid", "Kidney", "Liver", "Blood",
                  "Placenta", "Eye", "Heart", "Embryo", "Skeletal muscle",
                  "Brain"),
    "vis" = c("UMAP", "TSNE"),
    "ssclu" = c("ConstrainedKMeans", "SeededKMeans", "ConstrainedSeededKMeans"),
    "de" = c("TTest")
)
