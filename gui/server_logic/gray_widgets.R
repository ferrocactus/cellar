# Logic for graying out buttons when not applicable
gray_widgets <- function(input, output, session) {
    observeEvent(input$dim_method, {
        if(input$dim_method == "PCA" || input$dim_method == "Precomputed PCA"){
            shinyjs::enable(selector = "[type=radio][value=pca_auto]")
        } else {
            shinyjs::disable(selector = "[type=radio][value=pca_auto]")
            updateRadioButtons(
                session,
                "dim_options",
                selected = "pca_manual"
            )
        }
    })

    observeEvent(input$dim_options, {
        if(input$dim_options == "pca_manual"){
            shinyjs::enable('dim_n_components')
        } else {
            shinyjs::disable('dim_n_components')
        }
    })

    observeEvent(input$clu_method, {
        if(input$clu_method %in% options$clu_no_n_clusters){
            shinyjs::disable('clu_n_clusters')
            shinyjs::disable('eval_method')
        } else {
            shinyjs::enable('clu_n_clusters')
            shinyjs::enable('eval_method')
        }
    })
}
