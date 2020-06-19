align <- function(input, output, session, adata, selDatasetAlign,
                  replot, reset, relabel, resubset) {
    adataAlign <- reactiveVal(0)

    observeEvent(input$align_btn, {
        if (is_active(adata()) == FALSE) return()
        req(input$reference_dataset)

        withProgress(message = "Running SingleR", value = 0, {
            n <- 6
            incProgress(1 / n, detail = paste("Step: Reading data"))
            isolate(adataAlign(read_h5ad(input$reference_dataset$datapath)))
            #isolate(adataAlign(load_file(selDatasetAlign())))

            if (input$align_method == 'SingleR') {
                # transpose rows and cols for SingleR
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
            } else {
                msg <- cellar$safe(cellar$transfer_labels,
                    x = adata(),
                    ref = adataAlign(),
                    method = input$align_method,
                    inplace = TRUE
                )

                if (msg != 'good') {
                    showNotification(py_to_r(msg))
                    return()
                }
            }

            adataAlign(0)

            incProgress(1 / n, detail = "Converting names")
            msg <- cellar$safe(cellar$name_genes,
                x = adata(),
                inplace = TRUE
            )

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }

            incProgress(1 / n, detail = "Visualizing")
            msg <- cellar$safe(cellar$reduce_dim_vis,
                x = adata(),
                method = input$vis_method,
                dim = 2,
                use_emb = TRUE,
                inplace = TRUE,
                check_if_exists = TRUE)

            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }
        })

        replot(replot() + 1)
        reset(reset() + 1)
        relabel(relabel() + 1)
        resubset(resubset() + 1)
    })
}