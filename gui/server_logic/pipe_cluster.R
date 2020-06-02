# Needed by cluster_run
pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        eval_method, clu_n_clusters, vis_method) {
    withProgress(message = "Making plot", value = 0, {
        n <- 4
        incProgress(1 / n, detail = paste("Step: Reducing Data"))
        msg <- pipe()$run_step(step = 'dim', dim_method = dim_method,
                              dim_n_components = dim_n_components)
        if (msg != 'good') return(msg)
        incProgress(1 / n, detail = paste("Step: Clustering"))
        msg <- pipe()$run_step(step = 'clu', clu_method = clu_method,
                                eval_method = eval_method,
                                clu_n_clusters = clu_n_clusters)
        if (msg != 'good') return(msg)
        incProgress(1 / n, detail = paste("Step: Visualizing"))
        msg <- pipe()$run_step(step = 'vis', vis_method = vis_method)

        return(msg)
    })
}

pipe_sscluster <- function(pipe, new_labels, ssc_method, saved_clusters) {
    withProgress(message = "Running", value = 0, {
        incProgress(1 / 2, detail = paste("Constrained Clustering"))
        msg <- pipe()$run_step(step = "ssclu", ssclu_method = ssc_method,
                               ssclu_new_labels = new_labels,
                               saved_clusters = saved_clusters)
        return(msg)
    })
}

cluster_run <- function(input, output, session, pipe, selDataset, replot,
                        reset, newLabels) {
    # Clustering
    observeEvent(input$runconfigbtn, {
        if (pipe() == 0) {
            showNotification("Please load the data first.")
            return()
        }

        if (input$dim_options == "pca_auto")
            dim_n_components = 'knee'
        else
            dim_n_components = input$dim_n_components

        msg <- pipe_cluster(pipe, dim_method = input$dim_method,
                        dim_n_components = dim_n_components,
                        clu_method = input$clu_method,
                        eval_method = input$eval_method,
                        clu_n_clusters = input$clu_n_clusters,
                        vis_method = input$vis_method)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(1) # Notify that labels have changed
        reset(1) # notify changes
    })

    # Semi-supervised clustering
    observeEvent(input$ssclurun, {
        if (is.null(newLabels())) {
            showNotification("No labels have been updated")
            return()
        }

        msg <- pipe_sscluster(pipe, new_labels = newLabels(),
                              ssc_method = input$ssc_method,
                              saved_clusters = input$saved_clusters)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(1)
        reset(1)
    })
}
