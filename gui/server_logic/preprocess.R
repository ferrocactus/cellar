preprocess <- function(input, output, session, adata) {
    observeEvent(input$runpreprocessbtn, {
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        if (input$preprocess_method == 'Defaults') {
            withProgress(message = "Preprocessing", value = 0.5, {
                adata(cellar$preprocess(x = adata(), method='Scanpy (Defaults)'))
            })
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
