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
            n_components = input$n_components

        withProgress(message = "Please Wait", value = 0, {
            n <- 3
            incProgress(0 / n, detail = "Reducing Dimensionality")

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
                inplace = TRUE)

            incProgress(1 / n, detail = "Visualizing")
            cellar$reduce_dim_vis(
                x = adata(),
                method = input$vis_method,
                dim = 2,
                use_emb = TRUE,
                inplace = TRUE,
                check_if_exists = TRUE)
        })

        replot(replot() + 1) # Notify that labels have changed
        #reset(reset() + 1) # notify changes
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
        if (!py_has_attr(adata()$obs, 'labels')) return()

        output$cell_names_outp <- renderTable(width = "100%", {
            labels <- py_to_r(get_cluster_label_list(adata()))
            names <- py_to_r(get_cluster_name_list(adata()))
            df <- data.frame(as.character(labels), as.character(names))
            colnames(df) <- c("Cluster ID", "Name")
            isolate(relabel(0))
            return(df)
        })
    })
}