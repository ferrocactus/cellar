library(shiny)

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
                        accept = c(".json")
                    )
                )
            ),
            #div(
            #    class = "div_step div_download_cells",
            #    list(
            #        selectInput(
            #            ns("cell_subset_download"),
            #            "Select subset",
            #            choices=c("")
            #        ),
            #        downloadButton(
            #            ns("download_cells"),
            #            "Download Subset Data",
            #            class = "longbtn downloadbtn"
            #        )
            #    )
            #),
            div(
                class = "div_step div_download_plot",
                list(
                    splitLayout(
                        selectInput(
                            ns("plot_download_format"),
                            "Select format",
                            choices=c("SVG", "HTML", "JPEG", "PNG", "WEBP", "PDF")
                        ),
                        downloadButton(
                            ns("download_plot"),
                            "Download Plot",
                            class = "secondcol downloadbtn"
                        )
                    )
                )
            )
        )
    )
)}
