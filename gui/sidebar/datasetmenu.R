library(shiny)

datasetmenu <- menuItem(
    "Dataset",
    id = "datasetbtn",
    icon = icon("archive"),
    startExpanded = TRUE,
    menuSubItem(
        icon = NULL,
        list(
            radioButtons(
                "folder",
                "Select dataset group:",
                c(
                    "Uploaded Datasets" = "choose_temp",
                    "Server Datasets" = "choose_server"
                ),
                inline = TRUE
            ),

            conditionalPanel(
                condition = "input.folder == 'choose_temp'",
                selectInput(
                    "dataset",
                    "Choose uploaded dataset:",
                    choices = list.files("datasets")
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
            ),

            conditionalPanel(
                condition = "input.folder == 'choose_server'",
                selectInput(
                    "server_dataset",
                    "Choose server dataset:",
                    choices = list.files("server_datasets")
                )
            )
        )
    )
)
