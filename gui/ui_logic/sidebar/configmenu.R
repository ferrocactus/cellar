library(shiny)

source("gui/ui_logic/sidebar/options.R")

configmenu <- function(id, label="configmenu") {
ns = NS(id)
menuItem(
    "Clustering",
    id = "configbtn",
    icon = icon("cog"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            actionButton(
                ns("runconfigbtn"),
                "Run with current configuration",
                class="sidebtn longbtn"
            ),

            div(
                class = "div_step div_dim",
                list(
                    selectInput(
                        ns("dim_method"),
                        "Dimensionality reduction",
                        choices=options$dim
                    ),
                    splitLayout(
                        cellWidths = c("60%", "47%"),
                        radioButtons(
                            ns("dim_options"),
                            "Number of components",
                            choices = c(
                                "Automatic" = "pca_auto",
                                "Manual" = "pca_manual"
                            ),
                            inline = TRUE,
                            selected = 'pca_auto',
                        ),
                        div(
                            class = 'div_dim_n_components',
                            textInput(
                                ns("dim_n_components"),
                                label = NULL,
                                value = defaults$dim
                        ))
                    )
                ),
                HTML('<hr class="line">')
            ),

            div(
                class = "div_step div_clu",
                list(
                    selectInput(
                        ns("clu_method"),
                        "Clustering",
                        choices = options$clu
                    ),
                    splitLayout(
                        cellWidths = c("50%", "50%"),
                        textInput(
                            ns("clu_n_clusters"),
                            "Number of clusters",
                            value = defaults$clu
                        ),
                        selectInput(
                            ns("eval_method"),
                            "Evaluation:",
                            choices = options$eval
                        )
                    ),
                    conditionalPanel(
                        condition = "input.clu_method == 'Ensemble'",
                        ns = ns,
                        div(
                            class = "multicol",
                            checkboxGroupInput(
                                ns("ensemble_checkbox"),
                                label = NULL,
                                choices = c(
                                    options$clu_ensemble
                                )
                            )
                        )
                    )
                ),
                HTML('<hr class="line">')
            ),

            #div(
            #    class = "div_step div_mark",
            #    list(
            #        selectInput(
            #            "de_method",
            #            "DE",
            #            choices = options$de
            #        ),
            #        splitLayout(
            #            cellWidths = c("25%", "25%", "48%"),
            #            textInput(inputId = "mark_alpha",
            #                      label = "alpha",
            #                      value = defaults$mark_alpha),
            #            textInput(inputId = "mark_markers_n",
            #                      label = "markers",
            #                      value = defaults$mark_no),
            #            selectInput("mark_correction",
            #                        "Correction",
            #                        choices = options$correction)
            #        )
            #    ),
            #    HTML('<hr class="line">')
            #),

            #splitLayout(
            #    cellWidths = c("50%", "50%"),
            #    selectInput("con_convention",
            #                "Converter convention",
            #                choices = options$converter),
            #    textInput(inputId = "con_path",
            #              label = "Path (dont use)",
            #              value = '')
            #),

            #splitLayout(
            #    cellWidths = c("50%", "50%"),
            #    selectInput("ide_tissue",
            #                "Converter convention",
            #                choices = options$tissue),
            #    textInput(inputId = "ide_path",
            #              label = "Path (dont use)",
            #              value = "")
            #),

            div(
                class = "div_step div_vis",
                list(
                    selectInput(
                        ns("vis_method"),
                        "Visualization method:",
                        choices = options$vis
                    )
                ),
                HTML('<hr class="line">')
            ),

            div(
                class = "div_step div_ssclu",
                list(
                    actionButton(
                        ns("ssclurun"),
                        "Run constrained clustering",
                        class="sidebtn longbtn"
                    ),

                    splitLayout(
                        cellWidths = c("60%", "40%"),
                        selectInput(
                            ns("ssc_method"),
                            "Constrained clustering",
                            choices = options$ssclu
                        ),
                        textInput(
                            ns("savedClusters"),
                            label = "Preserve",
                            value = "1,2"
                        )
                    )
                )
            )
        )
    )
)}
