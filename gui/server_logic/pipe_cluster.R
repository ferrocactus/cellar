library(reticulate)

source_python("__init__.py")

# Needed by cluster_run
pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        eval_method, clu_n_clusters, vis_method) {
    withProgress(message = "Making plot", value = 0, {
        n <- 4
        incProgress(1/n, detail = paste("Step: Reducing Data"))
        pipe()$get_emb(method = dim_method, n_components = dim_n_components)
        incProgress(1/n, detail = paste("Step: Clustering"))
        pipe()$get_labels(method = clu_method, eval_method = eval_method,
                        n_clusters = clu_n_clusters)
        incProgress(1/n, detail = paste("Step: Visualizing"))
        pipe()$get_emb_2d(method = vis_method)
    })
}

cluster_run <- function(input, output, session, pipe, selDataset, setNames,
                        setPts, replot, remark, deButtons, deGenes, labelList) {
    observeEvent(input$runconfigbtn, {
        print(paste("Selected", selDataset()))

        if (pipe() == 0) {
            print("Initializing pipeline")
            pipe(Pipeline(x = selDataset()))
        } else if (pipe()$dataset != selDataset()) {
            print("Changing dataset")
            pipe()$restate(x = selDataset())
            setNames(c("None"))
            setPts(c(NA))
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

        replot(1) # Notify that labels have changed

        setNames(c("None"))
        setPts(c(NA))
        labelList(c())
        # Update sets
        for (i in 1:length(pipe()$n_clusters)) {
            #TODO replace cluster_i if rerun
            setNames(c(
                setNames(), (paste("Cluster_", as.character(i-1), sep = ""))))
            setPts(c(setPts(), list(which(pipe()$labels == (i-1)))))
            labelList(c(labelList(), i-1))
        }

        # Clear all analysis tabs
        if (length(deButtons()) > 0)
            for (i in 1:length(deButtons()))
                deButtons()[[i]]$destroy()
        deButtons(c())
        deGenes(c())
        output$DEbuttons = NULL
        output$genes = NULL
        output$GeneOntology = NULL
        output$KEGG = NULL
        output$Markers = NULL
        output$Msigdb = NULL
    })
}
