library(reticulate)

source_python("__init__.py")

load_dataset <- function(input, output, session, pipe, selDataset,
                         setNames, setPts) {
    observeEvent(input$load_dataset, {
        print(paste("Selected", selDataset()))

        withProgress(message = "Reading dataset", value = 0, {
            incProgress(1 / 2, detail = "Loading")
            if (pipe() == 0) {
                print("Initializing pipeline")
                pipe(Pipeline(x = selDataset()))
            } else if (pipe()$dataset != selDataset()) {
                print("Changing dataset")
                pipe()$restate(x = selDataset())
                setNames(c("None"))
                setPts(c(NA))
                output$Plot2 = NULL
            }
        })
    })
}
