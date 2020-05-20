library(shiny)

datasetmenu <- menuItem(
    "Dataset",
    id = "datasetbtn",
    icon = icon("archive"),
    startExpanded = TRUE,
    menuSubItem(
        icon = NULL,
        list(
            
            radioButtons("folder", "Select Which One to Use:",
                         c("Temporary Datasets" = "choose_temp",
                           "Server Datasets" = "choose_server"
                           )
                         ),
            
            selectInput(
                "dataset",
                "Choose from temporary datasets:",
                choices = list.files("datasets")
            ),
            
             selectInput(
                 "server_dataset",
                 "Choose from server datasets:",
                 choices = list.files("server_datasets")
             ),
            
            fileInput(
                "file1", "Choose CSV/h5ad File",
                multiple = FALSE,
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv", ".h5ad"
                )
            )
        )
    )
)
