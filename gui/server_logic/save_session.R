library(rjson)
library(tools)

save_session <- function(input, output, session, pipe, setNames, setPts,
                         deGenes, selDataset) {
    observe({
        output$download_sess <- downloadHandler(
            filename = function() {
                paste0(file_path_sans_ext(basename(selDataset())),
                       "_session.json")
            },
            content = function(file) {
                sess <- c()
                sess$'dataset' <- file_path_sans_ext(basename(selDataset()))
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
                    if (pipe()$has_emb()) {
                        sess$'x1_coordinates' <- pipe()$x_emb_2d[, 1]
                        sess$'x2_coordinates' <- pipe()$x_emb_2d[, 2]
                    }
                    if (pipe()$has_emb_2d())
                        sess$'labels' <- pipe()$labels
                }
                json_obj <- toJSON(sess)
                write(json_obj, file)
            }
        )
    })
}
