# Logic for graying out buttons when not applicable
misc <- function(input, output, session) {
    observeEvent(input$dim_method, {
        if(input$dim_method == "PCA" || input$dim_method == "Truncated SVD"){
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

    runjs("document.getElementById('ns-menu').onclick = function() {
       window.open('https://github.com/ferrocactus/cellar/blob/master/doc/cellar_guide.md', '_blank'); };")

    runjs("document.getElementById('ns-search').onclick = function() {
       window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + document.getElementById('ns-searchgene').value, '_blank'); };")
}

is_error <- function(msg, notify=FALSE, from_py=TRUE) {
    if (from_py) msg = py_to_r(msg)
    if (msg != 'good') {
        if (notify) showNotification(msg)
        return(TRUE)
    }
    return(FALSE)
}
