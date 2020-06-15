source_python("__init__.py")

pipe_align <- function(pipe, method, x_ref, col_ids_ref, labels_ref, vis_method,
                        key_maps) {
    withProgress(message = "Running Ingest", value = 0, {
        n <- 2
        incProgress(1 / n, detail = "Step: Annotating")
        msg <- pipe()$run_step('align', align_method = method, x_ref = x_ref,
                        col_ids_ref = col_ids_ref, labels_ref = labels_ref,
                        key_maps = key_maps)
        if (msg != 'good') return(msg)

        incProgress(1 / n, detail = paste("Step: Visualizing"))
        msg <- pipe()$run_step(step = 'vis', vis_method = vis_method)
        return(msg)
    })
}

align_run <- function(input, output, session, pipe, replot, relabel, reset,
                      selDatasetAlign, pipeAlign) {
    observeEvent(input$align_btn, {
        if (pipe() == 0) return()
        if (input$align_method == 'SingleR') return()

        req(input$upload_sess_align)

        pipeAlign(Pipeline(x = selDatasetAlign()))

        tryCatch({
            sess <- fromJSON(file = input$upload_sess_align$datapath)
            pipeAlign()$load_session(sess$'pipe_sess')

            msg <- pipe_align(pipe, input$align_method, pipeAlign()$x,
                            pipeAlign()$col_ids, pipeAlign()$labels,
                            vis_method = input$vis_method,
                            key_maps = pipeAlign()$key_maps)
        }, error = function(e) {
            print(e)
            showNotification("Error occurred in reading file.")
            return()
        })

        pipeAlign(0)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
    })
}
