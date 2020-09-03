scipy <- import('scipy')
anndata <- import('anndata')
library(densityClust)
source('gui/methods/epiConv/source/epiConv_functions.R')

cluster <- function(input, output, session, adata,
                    replot, reset, resubset) {
    observeEvent(input$runconfigbtn, {
        # Return if no adata loaded
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        # Determine n_components for dimensionality reduction
        if (input$dim_options == "pca_auto")
            n_components = 'knee'
        else {
            req(input$dim_n_components)
            n_components = input$dim_n_components
        }

        withProgress(message = "Please Wait", value = 0, {
            n <- 5

            if (py_to_r(is_sparse(adata())) && input$dim_method == 'ATAC') {
                if (!py_to_r(has_x_emb_sparse(adata(), input$dim_method, n_components))) {
                    atac_pipe(adata, n_components)
                }
            } else {
                incProgress(1 / n, detail = "Reducing Dimensionality")
                msg <- cellar$safe(cellar$reduce_dim,
                    x = adata(),
                    method = input$dim_method,
                    n_components = n_components,
                    inplace = TRUE,
                    check_if_exists = TRUE)

                if (is_error(msg)) return()

                incProgress(1 / n, detail = "Clustering")
                if (input$clu_method == 'Ensemble')
                    msg <- cellar$safe(cellar$cluster,
                        x = adata(),
                        method = input$clu_method,
                        eval_method = input$eval_method,
                        n_clusters = input$clu_n_clusters,
                        use_emb = TRUE,
                        inplace = TRUE,
                        check_if_exists = TRUE,
                        ensemble_methods = input$ensemble_checkbox)
                else if (input$clu_method == 'Leiden')
                    msg <- cellar$safe(cellar$cluster,
                        x = adata(),
                        method = input$clu_method,
                        use_emb = TRUE,
                        inplace = TRUE,
                        check_if_exists = FALSE,
                        resolution = input$leiden_resolution,
                        n_neighbors = input$leiden_neighbors)
                else
                    msg <- cellar$safe(cellar$cluster,
                        x = adata(),
                        method = input$clu_method,
                        eval_method = input$eval_method,
                        n_clusters = input$clu_n_clusters,
                        use_emb = TRUE,
                        inplace = TRUE,
                        check_if_exists = TRUE)

                if (is_error(msg)) return()

                incProgress(1 / n, detail = "Visualizing")
                msg <- cellar$safe(cellar$reduce_dim_vis,
                    x = adata(),
                    method = input$vis_method,
                    dim = 2,
                    use_emb = TRUE,
                    inplace = TRUE,
                    check_if_exists = TRUE)

                if (is_error(msg)) return()
            }

            incProgress(1 / n, detail = "Converting names")
            msg <- cellar$safe(cellar$name_genes,
                x = adata(),
                inplace = TRUE
            )

            if (is_error(msg)) return()
        })

        if (is_error(msg, notify = TRUE)) return()

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    # Semi-supervised clustering
    observeEvent(input$ssclurun, {
        req(adata())
        if (!py_has_attr(adata()$obs, 'labels')) return()

        withProgress(message = "Please Wait", value = 0, {
            n <- 2
            incProgress(1 / n, detail = "Clustering")
            msg <- cellar$safe(cellar$ss_cluster,
                x = adata(),
                method = input$ssc_method,
                use_emb = TRUE,
                inplace = TRUE,
                n_clusters = input$n_ss_clusters,
                preserved_labels = input$saved_clusters)
        })

        if (is_error(msg, notify = TRUE)) return()

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    # Merge clusters
    observeEvent(input$merge_clusters, {
        req(adata())
        req(input$clusters_to_merge)
        if (!py_has_attr(adata()$obs, 'labels')) return()

        msg <- cellar$safe(merge_clusters,
            adata = adata(),
            clusters = input$clusters_to_merge)

        if (is_error(msg, notify = TRUE)) return()

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })
}

atac_pipe <- function(adata, num.eigs) {
    n <- 5
    incProgress(1 / n, detail = "Preprocessing")

    mat <- py_to_r((adata()$T)$X)
    barcode <- as.character(py_to_r(adata()$obs_names$to_numpy()))
    rownames(mat) <- as.character(py_to_r(adata()$var_names$to_numpy()))
    colnames(mat) <- barcode
    lib_size <- lib.estimate(mat)
    freq <- freq.estimate(mat[,lib_size>1000])
    res_epiConv <- create.epiconv(
        meta.features = data.frame(barcode = barcode[lib_size>1000],
        lib_size = lib_size[lib_size>1000]))
    res_epiConv <- add.mat(obj = res_epiConv, x = mat[freq!=0, lib_size>1000], name = "peak")

    mat <- tfidf.norm(mat = res_epiConv@mat[["peak"]], lib_size = res_epiConv$lib_size)
    infv <- inf.estimate(mat[,sample(1:ncol(mat), size=500)], sample_size=0.125, nsim=30)

    adata(set_adata(mat, transpose=TRUE))

    sample_size <- floor(nrow(mat)/8)
    nsim <- 30
    sample_matrix <- lapply(1:nsim,function(x) sample(1:nrow(mat), size = sample_size))

    Smat <- matrix(0, res_epiConv@ncell, res_epiConv@ncell)
    for(i in sample_matrix){
        Smat <- Smat+epiConv.matrix(mat = mat[i,], inf_replace = infv)
    }
    Smat <- Smat/nsim

    incProgress(1 / n, detail = "Computing Similarity Matrix")
    res_epiConv <- add.similarity(obj = res_epiConv, x = Smat, name = "sampl_simp")

    Smat <- batch.blur(Smat = res_epiConv[["sampl_simp"]], batch = NULL, knn = 50)
    res_epiConv <- add.similarity(res_epiConv, x = Smat, name = "sampl_simp_denoise")

    incProgress(1 / n, detail = "Running UMAP")
    umap_settings <- umap::umap.defaults
    umap_settings$input <- "dist"
    umap_settings$n_components <- 2
    umap_res <- umap::umap(max(Smat)-Smat, config = umap_settings)$layout
    res_epiConv <- add.embedding(obj = res_epiConv, x = umap_res, name = "sampl_simp_denoise")

    incProgress(1 / n, detail = "Clustering")
    ncluster <- 10
    dclust_obj <- densityClust(res_epiConv@embedding[["sampl_simp_denoise"]],gaussian=T)
    rho_cut <- quantile(dclust_obj$rho,0.5)
    delta_cut <- sort(dclust_obj$delta[dclust_obj$rho>=rho_cut],decreasing=T)[ncluster+1]
    clust <- findClusters(dclust_obj,rho=rho_cut,delta=delta_cut)$clusters

    store_x_emb_d(adata(), umap_res, "UMAP")
    store_labels(adata(), clust, "densityClust")


    # if (num.eigs == 'knee') num.eigs = 50
    # x = scipy$sparse$csc_matrix(adata()$X)

    # # Random bins, we won't be needing them
    # bins <- GRanges(seqnames = "chr1",
    #               strand = c("+"),
    #               ranges = IRanges(start = c(1:dim(x)[2]), width = 3))

    # x.sp = createSnapFromBmat(x, barcodes=barcodes, bins=bins)
    # x.sp = runDiffusionMaps(x.sp, num.eigs=as.numeric(num.eigs))

    # plotDimReductPW(
    #   obj=x.sp,
    #   eigs.dims=1:50,
    #   point.size=0.3,
    #   point.color="grey",
    #   point.shape=19,
    #   point.alpha=0.6,
    #   down.sample=5000,
    #   pdf.file.name='EigenPlots.pdf',
    #   pdf.height=7,
    #   pdf.width=7
    # )

    # message(sprintf("Graph-based clustering\n"))
    # x.sp = runKNN(
    #   obj=x.sp,
    #   eigs.dims=1:20,
    #   k=15
    # )
    # x.sp=runCluster(
    #   obj=x.sp,
    #   tmp.folder='.',
    #   louvain.lib="R-igraph",
    #   seed.use=10
    # )

    # cellar$store_labels(adata(), x.sp@cluster, method='snapATAC')

    # x_emb = as.matrix(x.sp@smat@dmat)

    # return(x_emb)
}
