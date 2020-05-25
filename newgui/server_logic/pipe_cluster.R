library(reticulate)

source_python("__init__.py")
source("newgui/server_logic/selected_dataset.R")

pipe_cluster <- function(pipe, dim_method, dim_n_components, clu_method,
                        eval_method, clu_n_clusters, vis_method) {
    withProgress(message = "Making plot", value = 0, {
        n <- 4
        incProgress(1/n, detail = paste("Step: PCA"))
        pipe$get_emb(method = dim_method, n_components = dim_n_components)
        incProgress(1/n, detail = paste("Step: Clustering"))
        pipe$get_labels(method = clu_method, eval_method = eval_method,
                        n_clusters = clu_n_clusters)
        incProgress(1/n, detail = paste("Step: Visualizing"))
        pipe$get_emb_2d(method = vis_method)
    })
}

cluster_run <- function(input, output, session, selDataset, color) {
    observeEvent(input$runconfigbtn, {
        print(paste("Selected", selDataset()))

        if (TRUE) {
            print("Changing dataset")
            pipe <- Pipeline(x = selDataset())
        } else {
            # TODO
            print("Not changing dataset")
            pipe <- Pipeline(x = selDataset())
        }

        # Run clustering
        pipe_cluster(pipe, dim_method = input$dim_method,
                    dim_n_components = input$dim_n_components,
                    clu_method = input$clu_method,
                    eval_method = input$eval_method,
                    clu_n_clusters = input$clu_n_clusters,
                    vis_method = input$vis_method)

        # Update plot
        withProgress(message = "Making plot", value = 0, {
            incProgress(1/4, detail = paste("Step: Rendering plot"))
            if (input$color == "cluster") {
                plotcols = as.factor(pipe$labels)
            } else {
                plotcols = pipe$labels
            }

            output$plot <- renderPlotly({
                plot_ly(
                    x = pipe$x_emb_2d[,1], y = pipe$x_emb_2d[,2],
                    text = ~paste("Label: ", as.factor(pipe$labels)),
                    color = pipe$labels,
                    key = plotcols,
                    type = 'scatter',
                    mode = 'markers'
                ) %>% layout(dragmode = "lasso",
                            title = paste("Value of ", input$color, sep=""))
            })
        })

        return(pipe)
    })
}
