library(shiny)

downloadmenu <- menuItem(
    "Download",
    id = "downloadbtn",
    icon = icon("download"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_download_cells",
                list(
                    selectInput(
                        "cell_subset_download",
                        "Select subset",
                        choices=c("")
                    ),
                    downloadButton(
                        "download_cells",
                        "Download Subset Data",
                        class = "longbtn downloadbtn"
                    )
                )
            ),
            div(
                class = "div_step div_download_plot",
                list(
                    selectInput(
                        "plot_download_format",
                        "Select format",
                        choices=c("SVG", "PNG", "JPG")
                    ),
                    downloadButton(
                        "download_plot",
                        "Download Plot",
                        class = "longbtn downloadbtn"
                    )
                )
            )
        )
    )
)
