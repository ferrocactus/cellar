preprocess <- function(input, output, session, adata) {
    observeEvent(input$runpreprocessbtn, {
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        if (input$preprocess_method == 'Preprocess') {
            withProgress(message = "Preprocessing", value = 0.5, {
                filter_cells_min = as.numeric(input$filter_cells_min)
                filter_cells_max = as.numeric(input$filter_cells_max)
                filter_genes_min = as.numeric(input$filter_genes_min)
                filter_genes_max = as.numeric(input$filter_genes_max)
                normalize_total = as.numeric(input$normalize_total)
                if (input$exclude_highly_expressed == "Yes")
                    exclude_highly_expressed = TRUE
                else
                    exclude_highly_expressed = FALSE
                high_min_mean = as.numeric(input$high_min_mean)
                high_max_mean = as.numeric(input$high_max_mean)
                high_min_disp = as.numeric(input$high_min_disp)
                if (input$apply_log == "Yes")
                    apply_log = TRUE
                else
                    apply_log = FALSE
                scale = as.numeric(input$scale)

                filter_cells = list()
                if (filter_cells_min >= 0)
                    filter_cells[["run1"]] = list("min_genes" = filter_cells_min)
                if (filter_cells_max >= 0)
                    filter_cells[["run2"]] = list("max_genes" = filter_cells_max)

                filter_genes = list()
                if (filter_genes_min > -1)
                    filter_genes[["run1"]] = list("min_cells" = filter_genes_min)
                if (filter_genes_max > -1)
                    filter_genes[["run2"]] = list("max_cells" = filter_genes_max)

                normalize_total = list(
                    "target_sum" = normalize_total,
                    "exclude_highly_expressed" = exclude_highly_expressed
                )

                highly_variable_genes = list(
                    "min_mean" = high_min_mean,
                    "max_mean" = high_max_mean,
                    "min_disp" = high_min_disp
                )

                scale = list(
                    "max_value" = scale
                )

                a <- cellar$safe(cellar$preprocess,
                    x = adata(),
                    method='Scanpy',
                    filter_cells=filter_cells,
                    filter_genes=filter_genes,
                    normalize_total=normalize_total,
                    apply_log1p=apply_log,
                    highly_variable_genes=highly_variable_genes,
                    scale=scale)
            })

            if (py_to_r(is_str(a))) {
                showNotification(py_to_r(a))
                return()
            }

            adata(a)

        } else if (input$preprocess_method == 'sc-ATAC-seq') {
            withProgress(message = "Generating cell by gene matrix",
                         detail="This may take a while... (approx 3 mins)",
                         value = 0.4, {
                if (input$atac_strand == "No") {
                    stream_direction = FALSE
                    op_extend = NULL
                    max_op_extend = NULL
                } else {
                    stream_direction = TRUE
                    op_extend = input$atac_op_extend
                    max_op_extend = input$atac_max_op_extend
                }
                a <- cellar$safe(cellar$preprocess,
                    x = adata(),
                    method=input$preprocess_method,
                    operation=input$atac_operation,
                    extend=input$atac_extend,
                    max_extend=input$atac_max_extend,
                    stream_direction=stream_direction,
                    op_extend=op_extend,
                    max_op_extend=max_op_extend)
            })

            if (py_to_r(is_str(a))) {
                showNotification(py_to_r(a))
                return()
            }

            adata(a)
        }

        shape = py_to_r(adata()$shape)
        w = shape[[1]]
        h = shape[[2]]
        showNotification(paste0(
            "Finished preprocessing. New data has shape (", w, ", ", h, ")."))

    })
}
