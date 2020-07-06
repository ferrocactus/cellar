dataset <- function(input, output, session, adata, selDataset,
                    selDatasetAlign, fullreset) {
    ###########################################################################
    # Main Dataset
    toListen <- reactive({
        list(input$folder, input$uploaded_dataset, input$server_dataset)
    })

    # Stores selected dataset
    observeEvent(toListen(), {
        if (input$folder == 'user_uploaded')
            isolate(selDataset(as.character(input$uploaded_dataset)))
        else
            isolate(selDataset(as.character(input$server_dataset)))

        # Include path as well
        isolate(selDataset(
            paste0("datasets/", input$folder, "/", selDataset())))
    })

    ###########################################################################
    # Loads the dataset into an anndata object
    observeEvent(input$load_dataset, {
        print(paste("Selected", selDataset()))

        withProgress(message = "Please wait", value = 0, {
            incProgress(1 / 2, detail = "Reading dataset")
            isolate(adata(safe_load_file(selDataset())))
        })
        if (py_to_r(is_str(adata()))) {
            showNotification("Incorrect file format.")
            isolate(adata(0))
        } else {
            fullreset(fullreset() + 1)
            showNotification("Dataset loaded")
        }
    })

    ###########################################################################
    # Dataset used for alignment
    # toListenAlign <- reactive({
    #     list(input$folder_align, input$uploaded_dataset_align,
    #          input$server_dataset_align)
    # })

    # observeEvent(toListenAlign(), {
    #     if (input$folder_align == 'user_uploaded')
    #         isolate(selDatasetAlign(as.character(input$uploaded_dataset_align)))
    #     else
    #         isolate(selDatasetAlign(as.character(input$server_dataset_align)))

    #     # Include path as well
    #     isolate(selDatasetAlign(
    #         paste("datasets/", input$folder_align, "/",
    #                 selDatasetAlign(), sep = "")))
    # })
    ###########################################################################
}