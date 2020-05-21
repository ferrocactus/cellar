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
                    "Uploaded Datasets" = "user_uploaded",
                    "Server Datasets" = "hubmap"
                ),
                inline = TRUE
            ),

            conditionalPanel(
                condition = "input.folder == 'user_uploaded'",
                selectInput(
                    "uploaded_dataset",
                    "Choose dataset:",
                    choices = list.files("datasets/user_uploaded")
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
                condition = "input.folder == 'hubmap'",
                selectInput(
                    "hubmap_dataset",
                    "Choose dataset:",
                    choices = list.files("datasets/hubmap")
                )
            )
        )
    )
)
