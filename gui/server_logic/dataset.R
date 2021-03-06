
dataset <- function(input, output, session, adata, selDataset,
                    selDatasetAlign, fullreset, activeDataset) {
    ###########################################################################

    # udpate data tissue
    observeEvent(input$data_source, {
        updateSelectInput(
            session,
            "data_tissue",
            choices = list.files(paste0("datasets/server/", input$data_source))
        )
    })

    # update server data
    observeEvent(input$data_tissue, {
        updateSelectInput(
            session,
            "server_dataset",
            choices = list.files(paste0("datasets/server/", input$data_source, "/", input$data_tissue))
        )
    })

    # Main Dataset
    toListen <- reactive({
        list(input$folder, input$uploaded_dataset, input$server_dataset)
    })

    # Stores selected dataset
    observeEvent(toListen(), {
        if (input$folder == 'user_uploaded')
            isolate(selDataset(as.character(input$uploaded_dataset)))
        else
            isolate(selDataset(as.character(paste0(input$data_source, "/", input$data_tissue, "/", input$server_dataset))))

        # Include path as well
        isolate(selDataset(
            paste0("datasets/", input$folder, "/", selDataset())))
    })

    ###########################################################################
    # Loads the dataset into an anndata object
    observeEvent(input$load_dataset, {
        print(paste("Selected", selDataset()))
        activeDataset(tools::file_path_sans_ext(basename(selDataset())))

        withProgress(message = "Please wait", value = 0, {
            incProgress(1 / 2, detail = "Reading dataset")
            s=as.character(selDataset())
            if (substr(s,nchar(s)-3,nchar(s))=='_10x'){
                isolate(adata(read_10x(as.character(selDataset()))))

            }
            else{
                isolate(adata(safe_load_file(selDataset())))
            }
        })
        if (py_to_r(is_str(adata()))) {
            showNotification("Incorrect file format.")
            isolate(adata(NULL))
        } else {
            fullreset(fullreset() + 1)
            showNotification("Dataset loaded")
        }
    })

    ###########################################################################
    # Dataset used for alignment
    toListenAlign <- reactive({
        list(input$folder_align, input$uploaded_dataset_align,
             input$server_dataset_align)
    })

    observeEvent(toListenAlign(), {
        if (input$folder_align == 'user_uploaded')
            isolate(selDatasetAlign(as.character(input$uploaded_dataset_align)))
        else
            isolate(selDatasetAlign(as.character(input$server_dataset_align)))

        # Include path as well
        isolate(selDatasetAlign(
            paste("datasets/", input$folder_align, "/",
                    selDatasetAlign(), sep = "")))
    })
    ###########################################################################
}