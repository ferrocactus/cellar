runPipe <- function(pipe, input) {

    dim_method <- input$dim_method
    dim_n_components <- input$dim_n_components
    clu_method <- input$clu_method
    eval_method <- input$eval_method
    clu_n_clusters <- input$clu_n_clusters
    mark_method <- 'TTest'
    mark_alpha <- input$mark_alpha
    mark_markers_n <- input$mark_markers_n
    mark_correction <- input$mark_correction
    con_method <- 'Converter'
    con_convention <- input$con_convention
    #con_path <- input$con_path
    con_path <- 'markers/gene_id_name.csv'
    ide_method <- 'HyperGeom'
    #ide_path <- input$ide_path
    ide_path <- 'markers/cell_type_marker.json'
    ide_tissue <- input$ide_tissue
    vis_method <- input$vis_method
    ssc_method <- input$ssc_method

    if (clu_method == 'Scanpy') {
        df <- pipe$run_scanpy()
        labels <- df$labels()
        x_emb_2d <- df$x_emb_2d()
        df <- data.frame(x1 = x_emb_2d[, 1],
                         x2 = x_emb_2d[, 2],
                         y = labels)
    } else {
        withProgress(message = 'Making plot', value = 0, {
            # Number of times we'll go through the loop
            n <- 6

            incProgress(1/n, detail = paste("Step: PCA"))
            x_emb = pipe$get_emb(method=dim_method, n_components=dim_n_components)

            incProgress(1/n, detail = paste("Step: Clustering"))
            labels <- pipe$get_labels(method=clu_method,
                                    eval_method=eval_method,
                                    n_clusters=clu_n_clusters)

            incProgress(1/n, detail = paste("Step: Finding markers"))
            markers <- pipe$get_markers(method='TTest', alpha=mark_alpha,
                                        markers_n=mark_markers_n,
                                        correction=mark_correction)

            incProgress(1/n, detail = paste("Step: Converting names"))
            markers <- pipe$convert(method=con_method, convention=con_convention,
                                    path=con_path)

            incProgress(1/n, detail = paste("Step: Identifying cells"))
            markers <- pipe$identify(method="HyperGeom", path=ide_path,
                                    tissue=ide_tissue)

            incProgress(1/n, detail = paste("Step: Visualizing"))
            x_emb_2d <- pipe$get_emb_2d(x_emb, y=labels, method=vis_method)
            df <- data.frame(x1 = x_emb_2d[, 1],
                             x2 = x_emb_2d[, 2],
                             y = labels)
        })
    }
    return(df)
}

runSSClu <- function(pipe, labels, input) {

    dim_method <- input$dim_method
    dim_n_components <- input$dim_n_components
    eval_method <- input$eval_method
    mark_method <- 'TTest'
    mark_alpha <- input$mark_alpha
    mark_markers_n <- input$mark_markers_n
    mark_correction <- input$mark_correction
    con_method <- 'Converter'
    con_convention <- input$con_convention
    #con_path <- input$con_path
    con_path <- 'markers/gene_id_name.csv'
    ide_method <- 'HyperGeom'
    #ide_path <- input$ide_path
    ide_path <- 'markers/cell_type_marker.json'
    ide_tissue <- input$ide_tissue
    vis_method <- input$vis_method
    ssc_method <- input$ssc_method
    saved_clusters <- input$savedClusters
    saved_clusters <- as.numeric(strsplit(saved_clusters, ',')[[1]])

    withProgress(message = 'Please wait', value = 0, {
        # Number of times we'll go through the loop
        n <- 5

        #incProgress(1/n, detail = paste("Step: PCA"))
        #x_emb = pipe$get_emb(method=dim_method, n_components=dim_n_components)

        incProgress(1/n, detail = paste("Step: Clustering"))
        if (ssc_method == "SeededKMeans") {
            labels <- pipe$update(method=ssc_method,
                                  new_labels=labels)
        } else {
            labels <- pipe$update(method=ssc_method,
                                  new_labels=labels,
                                  saved_clusters=saved_clusters)
        }

        incProgress(1/n, detail = paste("Step: Finding markers"))
        markers <- pipe$get_markers(method='TTest', alpha=mark_alpha,
                                    markers_n=mark_markers_n,
                                    correction=mark_correction)

        incProgress(1/n, detail = paste("Step: Converting names"))
        markers <- pipe$convert(method=con_method, convention=con_convention,
                                path=con_path)

        incProgress(1/n, detail = paste("Step: Identifying cells"))
        markers <- pipe$identify(method="HyperGeom", path=ide_path,
                                 tissue=ide_tissue)

        incProgress(1/n, detail = paste("Step: Visualizing"))
        x_emb_2d <- pipe$get_emb_2d(y=labels, method=vis_method)
        df <- data.frame(x1 = x_emb_2d[, 1],
                         x2 = x_emb_2d[, 2],
                         y = labels)
    })
    return(df)
}
