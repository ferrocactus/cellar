history <- function(id, label = "history") {
    ns = NS(id)
    fluidRow(
        span(class = "hist",
            conditionalPanel(
                'true',
                ns = ns,
                dropdownButton(
                    tags$h4("Plot Options"),
                    splitLayout(
                        radioButtons(
                            ns("show_names"),
                            "Show cluster names in plot?",
                            c(
                                "No" = "dont_show_names",
                                "Yes" = "show_names"
                            ),
                            inline = TRUE
                        ),
                        radioButtons(
                            ns("theme_mode"),
                            "Select theme:",
                            c(
                                "Light Mode" = "light_mode",
                                "Dark Mode" = "dark_mode"
                            ),
                            inline = TRUE
                        )
                    ),
                    sliderInput(
                        ns("dot_size"),
                        "Select marker size",
                        min = 1, max = 30, value = 5
                    ),
                    sliderInput(
                        ns("plot_height"),
                        "Select plot height",
                        min = 400, max = 800, value = 600
                    ),
                    sliderInput(
                        ns("value_t"),
                        label="Gene expression thresholds",
                        min=0,max=10,
                        value=c(0, 10),
                        step=1
                    ),
                    conditionalPanel(
                        "input.color2 != 'None'",
                        ns = ns,
                        sliderInput(
                            ns("value_t_2"),
                            label="Gene 2 expression thresholds",
                            min=0,max=10,
                            value=c(0, 10),
                            step=1
                        )
                    ),
                    actionButton(
                        ns("gray_cells"),'Grayout cells out of range',
                        class = "longbtn-dropdown"
                    ),

                    selectInput(
                        ns("color_by"), "Color by key",
                        c("Clusters")
                    ),

                    size = 'sm',
                    circle = TRUE,
                    inputId = ns("appearance_dropdown"),
                    #status = "danger",
                    icon = icon("gear"),
                    width = "500px",
                    right = TRUE,
                    margin = '40px',
                    tooltip = tooltipOptions(title = "Configure Plot")
                ),
                dropdownButton(
                    tags$h4("Export plot or store in a new tab"),
                    splitLayout(
                        textInput(
                            ns("plot_download_scale"),
                            "Scale",
                            value=1
                        ),
                        textInput(
                            ns("plot_download_width"),
                            "Image width",
                            value=1000
                        ),
                        textInput(
                            ns("plot_download_height"),
                            "Image height",
                            value=1000
                        )
                    ),
                    splitLayout(
                        selectInput(
                            ns("plot_download_format"),
                            "Select format",
                            choices=c("PNG", "JPEG", "SVG", "HTML", "PDF", "WEBP")
                        ),
                        downloadButton(
                            ns("download_plot"),
                            "Download Plot",
                            class = "secondcol-dropdown"
                        )
                    ),
                    splitLayout(
                        actionButton(
                            ns("store_plot"), "Store plot in a new tab",
                            class = "store-plot"
                        ),
                        actionButton(
                            ns("delete_plot"), "Delete current plot tab",
                            class = "delete-plot"
                        )
                    ),
                    size = 'sm',
                    circle = TRUE,
                    inputId = ns("store_plot_dropdown"),
                    #status = "danger",
                    icon = icon("save"),
                    width = "400px",
                    right = TRUE,
                    margin = '40px',
                    tooltip = tooltipOptions(title = "Store Plot")
                )
            )
        )
    )
}
