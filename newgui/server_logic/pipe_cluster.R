library(reticulate)

source_python("__init__.py")
source("newgui/server_logic/selected_dataset.R")

pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        eval_method, clu_n_clusters) {
    withProgress(message = "Making plot", value = 0, {
        n <- 2
        incProgress(1/n, detail = paste("Step: PCA"))
        pipe$get_emb(method = dim_method, n_components = dim_n_components)
        incProgress(1/n, detail = paste("Step: Clustering"))
        pipe$get_labels(method = clu_method, eval_method = eval_method,
                        n_clusters = clu_n_clusters)
    })
}

cluster_run <- function(input, output, session, selDataset) {
    observeEvent(input$runconfigbtn, {
        print(paste("Selected", selDataset()))

        if (TRUE) {
            print("Changing dataset")
            pipe <- Pipeline(x = selDataset())
        } else {
            print("Not changing dataset")
            pipe <- Pipeline(x = selDataset())
        }

        pipe_cluster(pipe, dim_method = input$dim_method,
                    dim_n_components = input$dim_n_components,
                    clu_method = input$clu_method,
                    eval_method = input$eval_method,
                    clu_n_clusters = input$clu_n_clusters)
        return(pipe)
    })
}
