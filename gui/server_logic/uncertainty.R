

uncertain <- function(input, output, session, adata) {
    observeEvent(input$uncertain,{
        #showNotification("calculating uncertainty")
        
        withProgress(message = "Please Wait", value = 0, {
            incProgress(1 / 3, detail = "Data processing...")
            cellar$get_neighbors(adata())
            incProgress(1 / 3, detail = "Calculating uncertainty...")
            cellar$uncertainty(adata())
        })
        showNotification("finished")
    })
}

