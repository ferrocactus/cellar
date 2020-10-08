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
                         value = 0.5, {
                adata(cellar$preprocess(x = adata(), method=input$preprocess_method))
            })
        }
        shape = py_to_r(adata()$shape)
        w = shape[[1]]
        h = shape[[2]]
        showNotification(paste0(
            "Finished preprocessing. New data has shape (", w, ", ", h, ")."))

    })
}
