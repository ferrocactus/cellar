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
                            "Upload Dataset" = "user_uploaded",
                            "Annotated Datasets" = "server",
                            "Side Plot" = "side_plot"
                        ),
                        inline = FALSE
                    ),
                    conditionalPanel(
                        condition = "input.folder_align == 'user_uploaded'",
                        ns = ns,
                        fileInput(
                            ns("reference_dataset"),
                            "Upload Reference Dataset (h5ad)",
                            multiple = FALSE,
                            accept = c(".h5ad")
                        )%>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Choose the reference dataset."),
                                        placement= "bottom"
                                    )
                            )
                    ),
                    conditionalPanel(
                        condition = "input.folder_align == 'server'",
                        ns = ns,
                        selectInput(
                            ns("server_dataset_align"),
                            "Choose Reference Dataset:",
                            choices = list.files("datasets/annotated")
                        )
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
