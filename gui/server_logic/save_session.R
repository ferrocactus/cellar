save_session <- function(input, output, session, adata, replot,
                         remark, labelList, relabel, resubset,
                         fullreset, reinfo, activeDataset) {
    observe({
        output$download_sess <- downloadHandler(
            filename = function() {
                if (py_to_r(is_active(adata())) == FALSE) return()
                paste0(activeDataset(), "_cellar.h5ad")
            },
            content = function(path) {
                if (py_to_r(is_active(adata())) == FALSE) return()
                withProgress(message = "Saving File", value = 0, {
                    incProgress(1 / 2, detail = "Compressing")
                    #write_key(adata(), 'labelList', labelList())
                    write_h5ad(adata(), path, compression = 9)
                })
            }
        )
    })

    observeEvent(input$upload_sess, {
        req(input$upload_sess)
        activeDataset(tools::file_path_sans_ext(basename(input$upload_sess$datapath)))
        withProgress(message = "Loading Session", value = 0, {
            incProgress(1 / 3)
            fullreset(fullreset() + 1)
            updateSelectInput(
                session = session,
                inputId = "color",
                choices = c("Clusters"),
                selected = "Clusters")
            adata(read_h5ad(input$upload_sess$datapath))
            if (py_to_r(is_str(adata()))) {
                showNotification("Incorrect file format")
                isolate(adata(0))
                return()
            }

            incProgress(1 / 3)
            replot(replot() + 1)
            resubset(resubset() + 1)
            incProgress(1 / 3)
            remark(remark() + 1)
        })

        #labelList(py_to_r(adata()$uns[['labelList']]))
    })

    observeEvent(input$load_ann_dataset, {
        activeDataset(tools::file_path_sans_ext(input$annotated_datasets))
        datapath = paste0('datasets/annotated/', input$annotated_datasets)
        withProgress(message = "Loading Session", value = 0, {
            incProgress(1 / 3)
            fullreset(fullreset() + 1)
            updateSelectInput(
                session = session,
                inputId = "color",
                choices = c("Clusters"),
                selected = "Clusters")
            adata(read_h5ad(datapath))
            if (py_to_r(is_str(adata()))) {
                showNotification("Incorrect file format")
                isolate(adata(0))
                return()
            }

            incProgress(1 / 3)
            replot(replot() + 1)
            resubset(resubset() + 1)
            incProgress(1 / 3)
            remark(remark() + 1)
        })
    })
}
