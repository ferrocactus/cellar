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
                    "Server Datasets" = "server"
                ),
                inline = TRUE
            )

            ,
            conditionalPanel(
                condition = "input.folder == 'user_uploaded'",
                ns = ns,
                selectInput(
                    ns("uploaded_dataset"),
                    "Choose dataset:",
                    choices = c("default", list.files("datasets/user_uploaded"))
                )
                %>%
                    shinyInput_label_embed(
                        shiny::icon("info-circle") %>%
                            bs_embed_tooltip(
                                "Choose a specific dataset for analysis",
                                placement= "bottom"
                            )
                    )
                ,
                fileInput(
                    ns("file1"),
                    "Choose CSV/h5ad File",
                    multiple = FALSE,
                    accept = c(
                        "text/csv",
                        "text/comma-separated-values",
                        ".csv", ".h5ad"
                    )
                )
                %>%
                    shinyInput_label_embed(
                        shiny::icon("info-circle") %>%
                            bs_embed_tooltip(
                                "Upload a dataset",
                                placement= "bottom"
                            )
                    )
            ),

            conditionalPanel(
                condition = "input.folder == 'server'",
                ns = ns,
                selectInput(
                    ns("data_source"),
                    "Choose data center:",
                    choices = list.files("datasets/server")
                ),
                selectInput(
                    ns("data_tissue"),
                    "Choose tissue:",
                    choices = c("")
                ),
                selectInput(
                    ns("server_dataset"),
                    "Choose dataset:",
                    choices = c("")
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
