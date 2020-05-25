# Return the dataset that has been selected and store it in selDataset (reactive)
selected_dataset <- function(input, output, session) {
    selDataset <- reactiveVal(0)

    toListen <- reactive({
        list(input$folder, input$uploaded_dataset, input$hubmap_dataset)
    })

    observeEvent(toListen(), {
        if (input$folder == 'user_uploaded')
            selDataset(as.character(input$uploaded_dataset))
        else
            selDataset(as.character(input$hubmap_dataset))

        # Include path as well
        selDataset(paste(input$folder, "/", selDataset(), sep=""))
    })

    return(selDataset)
}
