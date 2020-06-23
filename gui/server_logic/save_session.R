save_session <- function(input, output, session, adata, replot,
                         remark, labelList, relabel, resubset,
                         fullreset) {
    observe({
        output$download_sess <- downloadHandler(
            filename = function() {
                if (is_active(adata()) == FALSE) return()
                paste0(adata()$uns[['dataset']], "_cellar.h5ad")
            },
            content = function(path) {
                if (is_active(adata()) == FALSE) return()
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
        adata(read_h5ad(input$upload_sess$datapath))
        if (py_to_r(is_str(adata()))) {
            showNotification("Incorrect file format")
            isolate(adata(0))
            return()
        }

        fullreset(fullreset() + 1)
        replot(replot() + 1)
        relabel(relabel() + 1)
        remark(remark() + 1)
        resubset(resubset() + 1)

        #labelList(py_to_r(adata()$uns[['labelList']]))
    })
}
