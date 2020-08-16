info <- function(input, output, session, adata, relabel, reinfo) {
    info_val <- reactiveValues()

    # Keep info tab active at all times
    output$cell_names_outp <- NULL
    output$clustering_info <- NULL
    outputOptions(output, "cell_names_outp", suspendWhenHidden = FALSE)
    outputOptions(output, "clustering_info", suspendWhenHidden = FALSE)
    # Start hidden
    shinyjs::hide("cell_names_outp")
    shinyjs::hide("clustering_info")

    # Show/Hide info
    observeEvent(input$collapse_cell_names, {
        shinyjs::toggle("cell_names_outp")
        shinyjs::toggle("clustering_info")
    })

    # Observe change of cluster names
    observe({
        if (relabel() < 1) return()
        isolate(relabel(0))

        req(adata())

        labels <- py_to_r(get_cluster_label_list(isolate(adata())))
        names <- py_to_r(get_cluster_name_list(isolate(adata())))
        df <- data.frame(as.character(labels), as.character(names))
        colnames(df) <- c("Cluster ID", "Label")

        tb <- df %>% addHtmlTableStyle(
            align='l', css.cell = "padding-right: 10em;") %>%
            htmlTable(caption = "Cluster Labels", rnames = FALSE)

        info_val$cellNames <- list(tb)
    })

    observe({
        if (reinfo() < 1) return()
        isolate(reinfo(0))

        req(adata())

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
                if (py_to_r(has_key_tri(isolate(adata()), 'uns', cat, tag)))
                    mx[i, 1] = paste(as.character(py_to_r(get_key_tri(
                        isolate(adata()), 'uns', cat, tag))), collapse=', ')
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
        info_val$configs <- list(tb)
    })

    return(info_val)
}