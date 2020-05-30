library(rjson)
library(tools)

save_session <- function(input, output, session, pipe, setNames, setPts,
                         deGenes, selDataset, plotObj, replot, remark) {
    observe({
        output$download_sess <- downloadHandler(
            filename = function() {
                paste0(file_path_sans_ext(basename(selDataset())),
                       "_session.json")
            },
            content = function(file) {
            withProgress( 1 / 2, detail = "Creating json file", value = 0, {
                sess <- c()
                pipe_sess <- pipe()$save_session()
                sess$'dataset' <- selDataset()
                sess$'dim_reduction_method' <- input$dim_method
                if (input$dim_options == 'pca_auto')
                    sess$'dim_n_components' <- 'Automatic'
                else
                    sess$'dim_n_components' <- input$dim_n_components
                sess$'clustering_method' <- input$clu_method
                sess$'clustering_evaluation' <- input$eval_method
                if (input$clu_method == 'Ensemble')
                    sess$'ensemble_methods' <- input$ensemble_checkbox
                sess$'visualization_method' <- input$vis_method
                sess$'de_genes_number' <- input$mark_markers_n
                sess$'ttest_alpha' <- input$alpha
                sess$'ttest_correction' <- input$correction
                if (length(setNames()) > 0)
                    sess$'set_names' <- setNames()
                if (length(setPts()) > 0)
                    sess$'set_points' <- setPts()
                if (length(deGenes()) > 0)
                    sess$'de_genes' <- deGenes()
                if (pipe() != 0) {
                    sess$'pipe_sess' <- pipe_sess
                }
                json_obj <- toJSON(sess)
                write(json_obj, file)
            })}
        )
    })

    observeEvent(input$upload_sess, {
        req(input$upload_sess)
        sess <- fromJSON(file = input$upload_sess$datapath)
        tryCatch({
            selDataset(sess$'dataset')
            updateSelectInput(session, 'dim_method',
                              selected = sess$'dim_reduction_method')
            if (sess$'dim_n_components' == 'Automatic') {
                updateRadioButtons(session, 'dim_options',
                                   selected = 'pca_auto')
            } else {
                updateRadioButtons(session, 'dim_options',
                                   selected = 'pca_manual')
                updateTextInput(session, 'dim_n_components',
                                value = sess$'dim_n_components')
            }
            updateSelectInput(session, 'clu_method',
                              selected = sess$'clustering_method')
            updateSelectInput(session, 'eval_method',
                              selected = sess$'clustering_evaluation')
            if (sess$'clustering_method' == 'Ensemble') {
                updateCheckboxGroupInput(session, 'ensemble_checkbox',
                              selected = sess$'ensemble_methods')
            }
            updateSelectInput(session, 'vis_method',
                              selected = sess$'visualization_method')
            updateSliderInput(session, 'mark_markers_n',
                            value = sess$'de_genes_number')
            updateTextInput(session, 'alpha',
                            value = sess$'ttest_alpha')
            updateTextInput(session, 'correction',
                            value = sess$'ttest_correction')
            if (!is.null(sess$'set_names')) {
                setNames(sess$'set_names')
            }
            if (!is.null(sess$'set_points')) {
                setPts(sess$'set_points')
            }
            if (!is.null(sess$'de_genes')) {
                deGenes(sess$'de_genes')
            }
            if (!is.null(sess$'pipe_sess')) {
                if (pipe() == 0) {
                    showNotification("Please load dataset first.")
                } else {
                    pipe()$load_session(sess$'pipe_sess')
                    replot(1)
                    if (pipe()$has('markers')) {
                        remark(1)
                    }
                }
            }
        }, error = function(e) {
            print(e)
            showNotification("Error occurred in reading file.")
        })
    })
}
