downloadmenu <- function(id, label="downloadmenu") {
    ns = NS(id)
    menuItem(
        "Export/Import",
        id = "downloadbtn",
        icon = icon("download"),
        startExpanded = FALSE,
        menuSubItem(
            icon = NULL,
            list(
                div(
                    class = "div_step div_download_session",
                    list(

                        downloadButton(
                            ns("download_sess"),
                            "Export Session",
                            class = "longbtn downloadbtn"
                        ),
                        # %>%
                        #     shinyInput_label_embed(
                        #         shiny_iconlink() %>%
                        #             bs_embed_tooltip( title = paste0("Download the current session for sharing or working later"), placement = "bottom")
                        #     ),  # doesn't work for download button
                        textOutput("Import Session"),
                        fileInput(
                            ns("upload_sess"),
                            "Import Session",
                            multiple = FALSE,
                            accept = c(".h5ad")
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip( title = paste0("Load a local session"), placement = "bottom")
                            ),
                        div(
                            class = "annotated_div",
                            list(
                                selectInput(
                                    ns("annotated_datasets"),
                                    "Choose Annotated Dataset",
                                    choices = c(list.files("datasets/annotated"))
                                ),
                                actionButton(
                                    ns("load_ann_dataset"),
                                    "Load Annotated Dataset",
                                    class = "sidebtn longbtn"
                                )
                            )
                        )
                    )
                ),
                div(
                    class = "div_step div_download_cells",
                    list(
                        splitLayout(
                            textInput(
                                ns("cell_subset_download"),
                                "Input Subset IDs",
                                #value = defaults$clu
                            ) %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip( title = paste0("Integer, list, or range specified by a dash, e.g., 1-7."), placement = "bottom")
                                ),
                            downloadButton(
                                ns("download_cells"),
                                "Download Subsets",
                                class = "secondcol downloadbtn",
                                value = 0
                                #style="width:50%"
                            )
                        )
                    )
                )
            )
        )
    )}
