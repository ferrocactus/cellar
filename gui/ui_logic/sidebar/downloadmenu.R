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
                        fileInput(
                            ns("upload_sess"),
                            "Import Session",
                            multiple = FALSE,
                            accept = c(".h5ad")
                        )
                    )
                ),
                div(
                    class = "div_step div_download_plot",
                    list(
                        splitLayout(
                            selectInput(
                                ns("plot_download_format"),
                                "Select format",
                                choices=c("PNG", "JPEG", "SVG", "HTML", "PDF", "WEBP")
                            ),
                            downloadButton(
                                ns("download_plot"),
                                "Download Plot",
                                class = "secondcol downloadbtn"
                            )
                        )
                    )
                ),
                div(
                    class = "div_step div_download_cells",
                    list(
                        textInput(
                            ns("cell_subset_download"),
                            "Input Subset IDs"
                            #value = defaults$clu
                        ),
                        downloadButton(
                            ns("download_cells"),
                            "Download Selected Subset",
                            class = "longbtn downloadbtn",
                            value = 0
                            #style="width:50%"
                        )
                    )
                )
            )
        )
    )}
