align <- function(input, output, session, adata, selDatasetAlign,
                  replot, reset, relabel, resubset, reinfo) {
    adataAlign <- reactiveVal(0)

    observeEvent(input$align_btn, {
        req(adata())
        if (input$folder_align == 'user_uploaded') {
            req(input$reference_dataset)
            isolate(selDatasetAlign(input$reference_dataset$datapath))
        } else {
            path = input$server_dataset_align
            path = paste0('datasets/annotated/', path)
            isolate(selDatasetAlign(path))
        }

        withProgress(message = "Running Label Transfer", value = 0, {
            n <- 5
            incProgress(1 / n, detail = paste("Step: Reading data"))
            isolate(adataAlign(safe_load_file(selDatasetAlign())))
            if (py_to_r(is_str(adataAlign()))) {
                showNotification("Incorrect file format.")
                isolate(adataAlign(0))
                return()
            }
            incProgress(1 / n, detail = paste("Step: Running Label Transfer"))
            if (input$align_method == 'SingleR') {
                # transpose rows and cols for SingleR
                x1 = t(py_to_r(get_x(adata())))
                rownames(x1) = py_to_r(get_var_names(adata()))
                colnames(x1) = py_to_r(get_obs_names(adata()))
                x2 = t(py_to_r(get_x(adataAlign())))
                rownames(x2) = py_to_r(get_var_names(adataAlign()))
                colnames(x2) = py_to_r(get_obs_names(adataAlign()))

                labels = py_to_r(get_labels(adataAlign()))
                if (labels == "No labels found") {
                    showNotification("No labels found. Please populate adata.obs['labels'] key.")
                    return()
                }

                print("Running SingleR.")
                pred <- SingleR(test = x1, ref = x2,
                                labels = labels)

                msg <- cellar$safe(store_labels,
                    adata = adata(),
                    labels = as.integer(pred$labels),
                    method = 'SingleR')

                if (msg != 'good') {
                    showNotification(py_to_r(msg))
                    adataAlign(0)
                    return()
                }

                msg <- cellar$safe(merge_cluster_names,
                    adata = adata(),
                    ref = adataAlign())
            } else {
                msg <- cellar$safe(cellar$transfer_labels,
                    x = adata(),
                    ref = adataAlign(),
                    method = input$align_method,
                    inplace = TRUE
                )
            }

            adataAlign(0)
            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }

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
        resubset(resubset() + 1)
    })
}