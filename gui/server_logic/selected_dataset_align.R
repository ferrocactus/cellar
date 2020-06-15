selected_dataset_align <- function(input, output, session, selDatasetAlign) {
    toListenAlign <- reactive({
        list(input$folder_align, input$uploaded_dataset_align,
             input$hubmap_dataset_align)
    })

    observeEvent(toListenAlign(), {
        if (input$folder_align == 'user_uploaded')
            selDatasetAlign(as.character(input$uploaded_dataset_align))
        else
            selDatasetAlign(as.character(input$hubmap_dataset_align))

        # Include path as well
        selDatasetAlign(paste(input$folder_align, "/", selDatasetAlign(), sep=""))
    })

}
