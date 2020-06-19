cluster <- function(input, output, session, adata, replot,
                    reset, relabel, resubset) {
    observeEvent(input$runconfigbtn, {
        # Return if no adata loaded
        if (is_active(adata()) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        # Determine n_components for dimensionality reduction
        if (input$dim_options == "pca_auto")
            n_components = 'knee'
        else
            n_components = input$dim_n_components

        withProgress(message = "Please Wait", value = 0, {
            n <- 5
            incProgress(1 / n, detail = "Reducing Dimensionality")
            cellar$reduce_dim(
                x = adata(),
                method = input$dim_method,
                n_components = n_components,
                inplace = TRUE,
                check_if_exists = TRUE)

            incProgress(1 / n, detail = "Clustering")
            cellar$cluster(
                x = adata(),
                method = input$clu_method,
                eval_method = input$eval_method,
                n_clusters = input$clu_n_clusters,
                use_emb = TRUE,
                inplace = TRUE,
                ensemble_methods = input$ensemble_checkbox)

            incProgress(1 / n, detail = "Visualizing")
            cellar$reduce_dim_vis(
                x = adata(),
                method = input$vis_method,
                dim = 2,
                use_emb = TRUE,
                inplace = TRUE,
                check_if_exists = TRUE)

            incProgress(1 / n, detail = "Converting names")
            cellar$name_genes(
                x = adata(),
                inplace = TRUE
            )
        })

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
        resubset(resubset() + 1)
    })

    # Semi-supervised clustering
    observeEvent(input$ssclurun, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        withProgress(message = "Please Wait", value = 0, {
            n <- 2
            incProgress(1 / n, detail = "Clustering")
            cellar$ss_cluster(
                x = adata(),
                method = input$ssc_method,
                use_emb = TRUE,
                inplace = TRUE,
                preserved_labels = input$saved_clusters)
        })

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
        resubset(resubset() + 1)
    })

    # Merge clusters
    observeEvent(input$merge_clusters, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()
        if (is.null(input$clusters_to_merge)) return()

        merge_clusters(adata(), input$clusters_to_merge)

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
        resubset(resubset() + 1)
    })

    # Show/Hide cluster names
    observeEvent(input$collapse_cell_names, {
        shinyjs::toggle("cell_names_outp")
    })

    # Observe change of cluster names
    observe({
        if (relabel() < 1) return()
        output$cell_names_outp <- renderTable(width = "100%", {
            isolate(relabel(0))
            if (has_key(adata(), 'uns', 'cluster_names') == FALSE) return()

            labels <- py_to_r(get_cluster_label_list(adata()))
            names <- py_to_r(get_cluster_name_list(adata()))
            df <- data.frame(as.character(labels), as.character(names))
            colnames(df) <- c("Cluster ID", "Name")
            return(df)
        })
    })
}