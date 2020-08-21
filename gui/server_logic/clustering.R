scipy <- import('scipy')

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
            incProgress(1 / n, detail = "Reducing Dimensionality")

            if (py_to_r(is_sparse(adata()))) {
                if (!py_to_r(has_x_emb_sparse(adata(), input$dim_method, n_components))) {
                    x_emb = diff_map_sparse(adata(), n_components)
                    store_x_emb(adata(), x_emb=x_emb, method=input$dim_method)
                }
            } else {
                msg <- cellar$safe(cellar$reduce_dim,
                    x = adata(),
                    method = input$dim_method,
                    n_components = n_components,
                    inplace = TRUE,
                    check_if_exists = TRUE)

                if (is_error(msg)) return()
            }

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

diff_map_sparse <- function(adata, num.eigs) {
    if (num.eigs == 'knee') num.eigs = 50
    x = scipy$sparse$csc_matrix(adata$X)
    barcodes = as.character(py_to_r(adata$obs_names$to_numpy()))
    # Random bins, we won't be needing them
    bins <- GRanges(seqnames = "chr1",
                  strand = c("+"),
                  ranges = IRanges(start = c(1:dim(x)[2]), width = 3))

    x.sp = createSnapFromBmat(x, barcodes=barcodes, bins=bins)
    x.sp = runDiffusionMaps(x.sp, num.eigs=as.numeric(num.eigs))
    x_emb = as.matrix(x.sp@smat@dmat)

    return(x_emb)
}