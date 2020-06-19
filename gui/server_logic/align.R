align <- function(input, output, session, adata, selDatasetAlign,
                  replot, reset, relabel, resubset) {
    adataAlign <- reactiveVal(0)

    observeEvent(input$align_btn, {
        if (is_active(adata()) == FALSE) return()

        isolate(adataAlign(load_file(selDatasetAlign())))

        if (input$align_method == 'SingleR') {
            # transpose rows and cols for SingleR
            withProgress(message = "Running SingleR", value = 0, {
                n <- 4
                incProgress(1 / n, detail = paste("Step: Reading data"))

                x1 = t(py_to_r(adata()$X))
                rownames(x1) = py_to_r(get_var_names(adata()))
                colnames(x1) = py_to_r(get_obs_names(adata()))
                x2 = t(py_to_r(adataAlign()$X))
                rownames(x2) = py_to_r(get_var_names(adataAlign()))
                colnames(x2) = py_to_r(get_obs_names(adataAlign()))

                incProgress(1 / n, detail = paste("Step: Annotating"))
                pred <- SingleR(test = x1, ref = x2,
                                labels = py_to_r(get_labels(adataAlign())))

                store_labels(adata(), as.integer(pred$labels), 'SingleR')

                incProgress(1 / n, detail = "Visualizing")
                cellar$reduce_dim_vis(
                    x = adata(),
                    method = input$vis_method,
                    dim = 2,
                    use_emb = TRUE,
                    inplace = TRUE,
                    check_if_exists = TRUE)

                adataAlign(0)
            })
        } else {
            # pipeAlign(Pipeline(x = selDatasetAlign()))

            # tryCatch({
            #     sess <- fromJSON(file = input$upload_sess_align$datapath)
            #     pipeAlign()$load_session(sess$'pipe_sess')

            #     msg <- pipe_align(pipe, input$align_method, pipeAlign()$x,
            #                     pipeAlign()$col_ids, pipeAlign()$labels,
            #                     vis_method = input$vis_method,
            #                     key_maps = pipeAlign()$key_maps)
            # }, error = function(e) {
            #     print(e)
            #     showNotification("Error occurred in reading file.")
            #     return()
            # })

            adataAlign(0)
        }

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
        resubset(resubset() + 1)
    })
}