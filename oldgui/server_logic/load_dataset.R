load_dataset <- function(input, output, session) {
    # Return selected dataset path
    if (input$folder == "user_uploaded") {
        dataset = as.character(input$uploadedDataset)
    } else if (input$folder == 'hubmap') {
        dataset = as.character(input$hubmapDataset)
    }

    return(paste(input$folder, "/", dataset, sep=""))
}
