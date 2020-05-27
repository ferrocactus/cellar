library(shiny)

gray_widgets <- function(input, output, session) {
    # Gray out automatic if PCA not selected
    observeEvent(input$dim_method, {
        if(input$dim_method == "PCA"){
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

    # Gray out n_components if manual not selected
    observeEvent(input$dim_options, {
        if(input$dim_options == "pca_manual"){
            shinyjs::enable('dim_n_components')
        } else {
            shinyjs::disable('dim_n_components')
        }
    })

    # Gray out eval method if clustering method does not require it
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
