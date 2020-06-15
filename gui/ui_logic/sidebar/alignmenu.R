alignmenu <- function(id, label="alignmenu") {
ns = NS(id)
menuItem(
    "Label Transfer",
    id = "alignbtn",
    icon = icon("arrows-alt"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_align",
                list(
                    selectInput(
                        ns("align_method"),
                        "Select label transfer method",
                        options$ali
                    ),
                    radioButtons(
                        ns("folder_align"),
                        "Select dataset group:",
                        c(
                            "Uploaded Datasets" = "user_uploaded",
                            "Server Datasets" = "hubmap"
                        ),
                        inline = TRUE
                    ),
                    conditionalPanel(
                        condition = "input.folder_align == 'user_uploaded'",
                        ns = ns,
                        selectInput(
                            ns("uploaded_dataset_align"),
                            "Choose dataset:",
                            choices = list.files("datasets/user_uploaded")
                        )
                    ),
                    conditionalPanel(
                        condition = "input.folder_align == 'hubmap'",
                        ns = ns,
                        selectInput(
                            ns("hubmap_dataset_align"),
                            "Choose dataset:",
                            choices = list.files("datasets/hubmap")
                        )
                    ),
                    fileInput(
                        ns("upload_sess_align"),
                        "Import Session",
                        multiple = FALSE,
                        accept = c(".json")
                    ),
                    actionButton(
                        ns("align_btn"),
                        "Run Label Transfer",
                        class = "longbtn sidebtn"
                    )
                )
            )
        )
    )
)
}
