
uncertain <- function(input, output, session, adata) {
    observeEvent(input$uncertain,{
        showNotification("calculating uncertainty")
        cellar$get_neighbors(adata())
        cellar$uncertainty(adata())
        showNotification("finished")
    })
}

