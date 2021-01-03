defaults <- list(
    "dim" = '40',
    "clu" = '(3, 5, 1)',
    "mark_alpha" = 0.05,
    "mark_no" = 200
)

options <- list(
    "dim" = c("PCA", "Incremental PCA", "Kernel PCA", "Truncated SVD",
              "Diffusion Map", "cisTopic", "MDS", "UMAP", "Isomap",
              "Spectral Embedding", "Feature Agglomeration"),
    "clu" = c("Leiden", "KMeans", "KMedoids", "Spectral", "Agglomerative",
              "GaussianMixture", "Ensemble", "Fixed"),
    "clu_ensemble" = c("All"="All", "KMeans"="KMeans", "KMedoids"="KMedoids",
                       "GaussianMixture"="GaussianMixture",
                       "Spectral"="Spectral", "Agglomerative"="Agglomerative",
                       "Leiden"="Leiden"),
    "clu_no_n_clusters" = c("DBSCAN", "Leiden"),
    "eval" = c("Silhouette", "DaviesBouldin", "CalinskiHarabasz"),
    "correction" = c("holm-sidak", "bonferroni", "sidak", "holm",
                     "simes-hochberg", "hommel", "fdr_bh", "fdr_by",
                     "fdr_tsbh", "fdr_tsbky", "None"),
    "converter" = c("id-to-name", "name-to-id"),
    "tissues" = c("all", "blood", "brain", "embryo", "eye", "heart", "kidney",
                  "large intestine", "liver", "lymph", "muscle", "other",
                  "placenta", "small intestine", "spleen", "stomach", "thymus",
                  "thyroid", "clusters", "user defined"),
    "vis" = c("UMAP", "TSNE", "Diffusion Map", "MDS", "PCA", "Kernel PCA",
              "Isomap", "Spectral Embedding", "Feature Agglomeration"),
    "ssclu" = c("ConstrainedKMeans", "SeededKMeans"),
    "de" = c("TTest", "Rank", "Wald"),
    "ali" = c("Scanpy Ingest", "SingleR"),
    "tissues" = c('Clusters','User defined')
)
