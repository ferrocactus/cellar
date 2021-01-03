cluster <- function(input, output, session, adata,
                    replot, reset, resubset) {
    observeEvent(input$reduce_dim_and_vis, {
        # Return if no adata loaded
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        withProgress(message = "Please Wait", value = 0, {
            n <- 4

            if (input$dim_method == 'cisTopic') {
                incProgress(1 / n, detail = "Running cisTopic. This may take a while...")

                if (input$dim_options == "pca_auto")
                    n_components = c(25, 30, 35, 40, 45, 50)
                else
                    n_components = input$dim_n_components

                x_emb = get_cell_by_topic(adata(), topic_list=n_components)
                store_x_emb(adata(), x_emb, method='cisTopic')
            } else {
                # Determine n_components for dimensionality reduction
                incProgress(1 / n, detail = "Reducing Dimensionality")

                if (input$dim_options == "pca_auto")
                    n_components = 'knee'
                else {
                    req(input$dim_n_components)
                    n_components = input$dim_n_components
                }

                msg <- cellar$safe(cellar$reduce_dim,
                    x = adata(),
                    method = input$dim_method,
                    n_components = n_components,
                    inplace = TRUE,
                    check_if_exists = TRUE)

                if (is_error(msg)) return()
            }

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

        if (!py_has_attr(adata()$obs, 'labels'))
            store_empty_labels(adata())

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    observeEvent(input$run_clustering, {
        # Return if no adata loaded
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }
        if (py_to_r(has_emb(adata())) == FALSE) {
            showNotification("Please run dimensionality reduction first.")
            return()
        }

        withProgress(message = "Please Wait", value = 0, {
            incProgress(1 / 2, detail = "Clustering")
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

get_cell_by_topic <- function(adata, topic_list) {
    # Used for scATAC-seq data
    # generate cell by topic matrix
    # bin names should be in  chr:1-1000 format
    # data matrix should be sparse
    sp <- import('scipy.sparse')
    topic_list = as.numeric(topic_list)

    a = adata$T # Transpose for cisTopic
    x = sp$csc_matrix(a$X)
    rownames(x) = as.character(py_to_r(correct_bin_names(a$obs_names$to_numpy())))
    colnames(x) = as.character(py_to_r(a$var_names$to_numpy()))

    cc = createcisTopicObject(x)
    print('Created cisTopic object.')
    print('Running LDA Models...')

    if (length(topic_list) >= 3)
        nCores = min(8, length(topic_list))
    else
        nCores = length(topic_list)

    cc = runWarpLDAModels(
        cc,
        topic=topic_list,
        nCores=min(4, length(topic_list)),
        iterations=150,
        addModels=FALSE)

    if (length(topic_list) < 3)
        cc = selectModel(cc, type="maximum")
    else # type=derivative
        cc = selectModel(cc, keepModels=FALSE)

    # Normalized topic by cell matrix
    cellassign <- modelMatSelection(cc, 'cell', 'Probability')
    return(t(cellassign))
}