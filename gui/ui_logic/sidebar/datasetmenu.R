datasetmenu <- function(id, label='datasetmenu') {
ns = NS(id)
menuItem(
    "Dataset",
    id = "datasetbtn",
    icon = icon("archive"),
    startExpanded = TRUE,
    menuSubItem(
        icon = NULL,
        list(
            radioButtons(
                ns("folder"),
                "Select dataset group:",
                c(
                    "Uploaded Datasets" = "user_uploaded",
                    "Server Datasets" = "hubmap"
                ),
                inline = TRUE
            ),

            conditionalPanel(
                condition = "input.folder == 'user_uploaded'",
                ns = ns,
                selectInput(
                    ns("uploaded_dataset"),
                    "Choose dataset:",
                    choices = list.files("datasets/user_uploaded")
                ),
                fileInput(
                    ns("file1"),
                    "Choose CSV/h5ad File",
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
                ns = ns,
                selectInput(
                    ns("hubmap_dataset"),
                    "Choose dataset:",
                    choices = list.files("datasets/hubmap")
                )
            ),

            actionButton(
                ns("load_dataset"),
                "Load Dataset",
                class = "sidebtn longbtn"
            )
        )
    )
)}
