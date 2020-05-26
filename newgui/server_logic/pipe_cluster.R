library(reticulate)

source_python("__init__.py")

pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        eval_method, clu_n_clusters, vis_method) {
    withProgress(message = "Making plot", value = 0, {
        n <- 4
        incProgress(1/n, detail = paste("Step: Reducing Dimensionality"))
        pipe()$get_emb(method = dim_method, n_components = dim_n_components)
        incProgress(1/n, detail = paste("Step: Clustering"))
        pipe()$get_labels(method = clu_method, eval_method = eval_method,
                        n_clusters = clu_n_clusters)
        incProgress(1/n, detail = paste("Step: Visualizing"))
        pipe()$get_emb_2d(method = vis_method)
    })
}

cluster_run <- function(input, output, session, pipe, selDataset,
                        setNames, setPts, replot) {
    observeEvent(input$runconfigbtn, {
        print(paste("Selected", selDataset()))

        if (pipe() == 0) {
            print("Initializing pipeline")
            pipe(Pipeline(x = selDataset()))
        } else if (pipe()$dataset != selDataset()) {
            print("Changing dataset")
            pipe()$restate(x = selDataset())
        } else {
            print("Not changing dataset")
        }

        # Run clustering
        pipe_cluster(pipe, dim_method = input$dim_method,
                    dim_n_components = input$dim_n_components,
                    clu_method = input$clu_method,
                    eval_method = input$eval_method,
                    clu_n_clusters = input$clu_n_clusters,
                    vis_method = input$vis_method)

        # Update plot
        replot(replot() + 1)

        # Update sets
        for (i in 1:length(pipe()$n_clusters)) {
            #TODO replace cluster_i if rerun
            setNames(c(
                setNames(), (paste("Cluster_", as.character(i-1), sep = ""))))
            setPts(c(setPts(), list(which(pipe()$labels == (i-1)))))
        }
        # TODO
        # clear analysis tabs
    })
}
