cluster <- function(input, output, session, adata, replot,
                    reset, relabel, resubset, reinfo,
                    cellNamesTb, infoTb) {
    # Keep info tab active at all times
    output$cell_names_outp <- NULL
    output$clustering_info <- NULL
    outputOptions(output, "cell_names_outp", suspendWhenHidden = FALSE)
    outputOptions(output, "clustering_info", suspendWhenHidden = FALSE)
    # Start hidden
    shinyjs::hide("cell_names_outp")
    shinyjs::hide("clustering_info")

    observeEvent(input$runconfigbtn, {
        # Return if no adata loaded
        if (py_to_r(is_active(adata())) == FALSE) {
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
            msg <- cellar$safe(cellar$reduce_dim,
                x = adata(),
                method = input$dim_method,
                n_components = n_components,
                inplace = TRUE,
                check_if_exists = TRUE)

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
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
            else
                msg <- cellar$safe(cellar$cluster,
                    x = adata(),
                    method = input$clu_method,
                    eval_method = input$eval_method,
                    n_clusters = input$clu_n_clusters,
                    use_emb = TRUE,
                    inplace = TRUE,
                    check_if_exists = TRUE)

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }

            incProgress(1 / n, detail = "Visualizing")
            msg <- cellar$safe(cellar$reduce_dim_vis,
                x = adata(),
                method = input$vis_method,
                dim = 2,
                use_emb = TRUE,
                inplace = TRUE,
                check_if_exists = TRUE)

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }

            incProgress(1 / n, detail = "Converting names")
            msg <- cellar$safe(cellar$name_genes,
                x = adata(),
                inplace = TRUE
            )

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }
        })

        if (msg != 'good') {
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    # Semi-supervised clustering
    observeEvent(input$ssclurun, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()
        req(input$saved_clusters)
        req(input$n_ss_clusters)

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

        if (msg != 'good') {
            showNotification(py_to_r(msg))
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    # Merge clusters
    observeEvent(input$merge_clusters, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()
        req(input$clusters_to_merge)

        msg <- cellar$safe(merge_clusters,
            adata = adata(),
            clusters = input$clusters_to_merge)

        if (msg != 'good') {
            showNotification(py_to_r(msg))
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })

    # Show/Hide cluster names
    observeEvent(input$collapse_cell_names, {
        shinyjs::toggle("cell_names_outp")
        shinyjs::toggle("clustering_info")
    })

    # Observe change of cluster names
    observe({
        if (relabel() < 1) return()

        labels <- py_to_r(get_cluster_label_list(adata()))
        names <- py_to_r(get_cluster_name_list(adata()))
        df <- data.frame(as.character(labels), as.character(names))
        colnames(df) <- c("Cluster ID", "Label")

        tb <- df %>% addHtmlTableStyle(
            align='l', css.cell = "padding-right: 10em;") %>%
            htmlTable(caption = "Cluster Labels", rnames = FALSE)
        isolate(cellNamesTb(c(cellNamesTb(), list(tb))))
        isolate(relabel(0))
    })

    observe({
        if (reinfo() < 1) return()

        mx <- matrix(ncol = 1, nrow = 8)
        i <- 1
        cats = list(
            "dim_reduction_info" = c(
                "method", "n_components", "n_components_used"),
            "cluster_info" = c(
                "method", "n_clusters", "n_clusters_used", "eval_method"),
            "visualization_info_2d" = c(
                "method"))

        for (cat in names(cats)) {
            for (tag in get(cat, cats)) {
                if (py_to_r(has_key_tri(adata(), 'uns', cat, tag)))
                    mx[i, 1] = paste(as.character(py_to_r(get_key_tri(
                        adata(), 'uns', cat, tag))), collapse=', ')
                else
                    mx[i, 1] = "NA"
                i = i+1
            }
        }

        rownames(mx) <- c(
            "Method",
            "Number of Components",
            "Number of Components Used",
            "Method",
            "Number of Clusters",
            "Clusters Used",
            "Evaluation Method",
            "Method"
        )

        colnames(mx) <- c("Value")

        tb <- mx %>% addHtmlTableStyle(
            align='l',
            css.cell = "padding-right: 10em;"
            ) %>%
            htmlTable(
                caption = "General Info",
                rgroup = c("Dimensionality Reduction", "Clustering", "Visualization"),
                n.rgroup = c(3, 4, 1)
            )
        isolate(infoTb(c(infoTb(), list(tb))))
            isolate(reinfo(0))
    })
}