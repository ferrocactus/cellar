library(shiny)

source_python('src/utils/read.py') # upload_file

nameDialogInput <- function(id, code = 0) {
    # Show dialog box when user inputs dataset name
    ns <- NS(id)
    modalDialog(
        textInput(ns("datasetName"), "Enter dataset name."),
        if (code == 1) {
            div(tags$b("Dataset exists.", style = "color: red;"))
        } else if (code == 2) {
            div(tags$b("Name should not be empty.", style = "color: red;"))
        },
        footer = tagList(actionButton(ns("ok"), "OK"))
    )
}

datasetExists <- function(datasetName, path) {
    # Check if dataset exists in path
    files <- list.files(path)

    if (length(files) > 0)
        for (i in 1:length(files))
            if (tools::file_path_sans_ext(files[i]) == datasetName)
                return(TRUE)
    return(FALSE)
}

writeDataset <- function(datasetName, inputPath) {
    # File upload
    withProgress(message = 'Please wait...', value = 0, {
        incProgress(1/2, detail = paste("Processing File"))
        upload_file(input_path = inputPath, output_name = datasetName)
        incProgress(2/2, detail = paste("Saving File"))
    })
}

upload <- function(input, output, session) {
    # Prompt user for dataset name. The name should not be empty
    # and there should be no dataset with the same name.
    observeEvent(input$uploadedFile, {
        req(input$uploadedFile)
        showModal(nameDialogInput(id = "ns_upload", code = 0))
    })

    observeEvent(input$ok, {
        # Check if dataset exists
        datasetName <- input$datasetName

        if (datasetName == "") {
            showModal(nameDialogInput(id = "ns_upload", code = 2))
        } else if (datasetExists(datasetName, path = "datasets/user_uploaded")) {
            showModal(nameDialogInput(id = "ns_upload", code = 1))
        } else {
            removeModal()
            writeDataset(datasetName = datasetName,
                        inputPath = input$uploadedFile$datapath)
            updateSelectInput(
                session = session,
                inputId = "uploadedDataset",
                label = "Choose a dataset:",
                choices = list.files("datasets/user_uploaded"),
                showNotification("Dataset uploaded")
            )
        }
    })
}
