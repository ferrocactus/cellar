align <- function(input, output, session, adata, selDatasetAlign,
                  replot, reset, relabel, resubset, reinfo,
                  second_plot_path, double_plot) {
    adataAlign <- reactiveVal(0)

    observeEvent(input$align_btn, {
        if (py_to_r(is_active(adata())) == FALSE) {
            showNotification("Please load the data first.")
            return()
        }

        if (input$folder_align == 'user_uploaded') {
            req(input$reference_dataset)
            isolate(selDatasetAlign(input$reference_dataset$datapath))
        } else if (input$folder_align == 'server') {
            path = input$server_dataset_align
            path = paste0('datasets/annotated/', path)
            isolate(selDatasetAlign(path))
        } else if (input$folder_align == 'side_plot') {
            if (double_plot() == FALSE) {
                showNotification("No Side Plot found.")
                return()
            }
            path = isolate(second_plot_path())
            isolate(selDatasetAlign(path))
        } else {
            showNotification("Dataset Group not found.")
            return()
        }

        withProgress(message = "Running Label Transfer", value = 0, {
            n <- 5
            if (input$align_method == 'SingleR') {
                incProgress(1 / n, detail = paste("Step: Reading data"))
                isolate(adataAlign(safe_load_file(selDatasetAlign())))
                if (py_to_r(is_str(adataAlign()))) {
                    msg <- "Incorrect file format."
                    return()
                }

                incProgress(1 / n, detail = paste("Step: Running SingleR"))
                # transpose rows and cols for SingleR
                labels = py_to_r(get_labels(adataAlign()))
                if (labels == "No labels found") {
                    msg <- "No labels found. Please populate adata.obs['labels'] key."
                    return()
                }

                x1 = t(py_to_r(get_x(adata())))
                rownames(x1) = py_to_r(get_var_names(adata()))
                colnames(x1) = py_to_r(get_obs_names(adata()))
                x2 = t(py_to_r(get_x(adataAlign())))
                rownames(x2) = py_to_r(get_var_names(adataAlign()))
                colnames(x2) = py_to_r(get_obs_names(adataAlign()))

                print("Running SingleR.")
                pred <- SingleR(test = x1, ref = x2,
                                labels = labels)

                msg <- cellar$safe(store_labels,
                    adata = adata(),
                    labels = as.integer(pred$labels),
                    method = 'SingleR')

                if (!is_error(msg)) {
                    msg <- cellar$safe(merge_cluster_names,
                        adata = adata(),
                        ref = adataAlign())
                }
            } else {
                incProgress(1 / n, detail = paste("Step: Running Ingest. This may take a while."))
                msg <- cellar$safe(cellar$transfer_labels,
                    x = adata(),
                    ref = selDatasetAlign(),
                    method = input$align_method,
                    inplace = TRUE
                )
            }

            isolate(adataAlign(0))
            if (is_error(msg)) return()

            incProgress(1 / n, detail = "Converting names")
            msg <- cellar$safe(cellar$name_genes,
                x = adata(),
                inplace = TRUE
            )

            if (is_error(msg)) return()

            incProgress(1 / n, detail = "Visualizing")
            msg <- cellar$safe(cellar$reduce_dim_vis,
                x = adata(),
                method = input$vis_method,
                dim = 2,
                use_emb = TRUE,
                inplace = TRUE,
                check_if_exists = TRUE)

            if (is_error(msg)) return()
        })

        if (is_error(msg, notify=TRUE)) return()

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
    })
}