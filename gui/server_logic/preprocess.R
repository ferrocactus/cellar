preprocess <- function(input, output, session, adata) {
    observeEvent(input$runpreprocessbtn, {
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        if (input$preprocess_method == 'Defaults') {
            withProgress(message = "Preprocessing", value = 0.5, {
                adata(cellar$preprocess(x = adata()))
            })
            shape = py_to_r(adata()$shape)
            w = shape[[1]]
            h = shape[[2]]
            showNotification(paste0(
                "Finished preprocessing. New data has shape (", w, ", ", h, ")."))
        }
    })
}