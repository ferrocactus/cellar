# Needed by cluster_run
pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        ensemble_methods, eval_method, clu_n_clusters,
                        vis_method) {
    withProgress(message = "Making plot", value = 0, {
        n <- 4
        incProgress(1 / n, detail = paste("Step: Reducing Data"))
        msg <- pipe()$run_step(step = 'dim', dim_method = dim_method,
                              dim_n_components = dim_n_components)
        if (msg != 'good') return(msg)
        incProgress(1 / n, detail = paste("Step: Clustering"))
        msg <- pipe()$run_step(step = 'clu', clu_method = clu_method,
                                eval_method = eval_method,
                                clu_n_clusters = clu_n_clusters,
                                ensemble_methods = ensemble_methods)
        if (msg != 'good') return(msg)
        incProgress(1 / n, detail = paste("Step: Visualizing"))
        msg <- pipe()$run_step(step = 'vis', vis_method = vis_method)

        return(msg)
    })
}

pipe_sscluster <- function(pipe, ssc_method, saved_clusters) {
    withProgress(message = "Running", value = 0, {
        incProgress(1 / 2, detail = paste("Constrained Clustering"))
        msg <- pipe()$run_step(step = "ssclu", ssclu_method = ssc_method,
                               saved_clusters = saved_clusters)
        return(msg)
    })
}

cluster_run <- function(input, output, session, pipe, selDataset, replot,
                        reset, relabel) {
    # Semi-supervised clustering
    observeEvent(input$ssclurun, {
        if (pipe() == 0) return()
        if (!pipe()$has('labels')) return()
        msg <- pipe_sscluster(pipe, ssc_method = input$ssc_method,
                              saved_clusters = input$saved_clusters)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
    })

    # Merge clusters
    observeEvent(input$merge_clusters, {
        if (pipe() == 0) return()
        if (!pipe()$has('labels')) return()
        if (is.null(input$clusters_to_merge)) return()

        msg <- pipe()$run_step(step = 'merge_clusters',
                            clusters = input$clusters_to_merge)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
    })
}
