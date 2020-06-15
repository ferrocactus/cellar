singler <- function(input, output, session, pipe, replot, relabel, reset,
                    selDatasetAlign, pipeAlign) {
    observeEvent(input$align_btn, {
        if (pipe() == 0) return()
        if (input$align_method != 'SingleR') return()

        req(input$upload_sess_align)

        # transpose rows and cols for SingleR
        withProgress(message = "Running SingleR", value = 0, {
            n <- 3
            incProgress(1 / n, detail = paste("Step: Reading data"))

            pipeAlign(Pipeline(x = selDatasetAlign()))

            tryCatch({
                sess <- fromJSON(file = input$upload_sess_align$datapath)
                pipeAlign()$load_session(sess$'pipe_sess')
            }, error = function(e) {
                print(e)
                showNotification("Error occurred in reading file.")
                return()
            })

            x1 = t(pipe()$x)
            rownames(x1) = pipe()$col_ids
            colnames(x1) = pipe()$row_ids
            x2 = t(pipeAlign()$x)
            rownames(x2) = pipeAlign()$col_ids
            colnames(x2) = pipeAlign()$row_ids

            incProgress(1 / n, detail = paste("Step: Annotating"))
            pred <- SingleR(test = x1, ref = x2, labels = pipeAlign()$labels)
            pipe()$set_labels(pred$labels)

            incProgress(1 / n, detail = paste("Step: Visualizing"))
            msg <- pipe()$run_step(step = 'vis', vis_method = input$vis_method)

            pipeAlign(0)

            if (msg != 'good') {
                showNotification(msg)
                return()
            }
        })

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
    })
}
