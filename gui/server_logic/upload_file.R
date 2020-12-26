# File upload logic
source_python("gui/server_logic/write_file.py")

datasetExists <- function(dataset, path) {
    # Check if dataset exists in path
    files <- list.files(path)

    if (length(files) > 0)
        for (i in 1:length(files))
            if (tools::file_path_sans_ext(files[i]) == dataset)
                return(TRUE)
    return(FALSE)
}

writeDataset <- function(path, name) {
    # File upload
    withProgress(
        message = 'Please wait...',
        value = 0, {
            incProgress(1/2, detail = paste("Processing file"))
            write_file_py(name, path)
            incProgress(2/2, detail = paste("Saving file"))
        }
    )
}

upload_file <- function(input, output, session) {
    ns <- session$ns
    observeEvent(input$file1, {
        req(input$file1)
        showModal(modalDialog(
            textInput(ns("dataset_fname"), "Enter dataset name."),
            footer = tagList(actionButton(ns("okdataset"), "OK")),
            fade = FALSE
        ))
    })

    observeEvent(input$okdataset, {
        # Check if dataset exists
        fname <- input$dataset_fname

        if (fname == "") {
            showModal(modalDialog(
                textInput(ns("dataset_fname"), "Enter dataset name."),
                div(tags$b("Name should not be empty.", style = "color: red;")),
                footer = tagList(actionButton(ns("okdataset"), "OK")),
                fade = FALSE
            ))
        } else if (datasetExists(fname, path="datasets/user_uploaded")) {
            showModal(modalDialog(
                textInput(ns("dataset_fname"), "Enter dataset name."),
                div(tags$b("Dataset exists.", style = "color: red;")),
                footer = tagList(actionButton(ns("okdataset"), "OK")),
                fade = FALSE
            ))
        } else {
            removeModal()
            print(paste("Writing dataset", fname))
            writeDataset(input$file1$datapath, fname)
            updateSelectInput(
                session = session,
                inputId = "uploaded_dataset",
                label = "Choose a dataset:",
                choices = c("default", list.files("datasets/user_uploaded")),
            )
            showNotification("Dataset uploaded")
        }
        
        
    })
}
