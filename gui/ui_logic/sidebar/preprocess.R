preprocessmenu <- function(id, label="preprocessmenu") {
ns = NS(id)
menuItem(
    "Preprocessing",
    id = "preprocessingbtn",
    icon = icon("filter"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_preprocess",
                list(
                    selectInput(
                        ns("preprocess_method"),
                        "Preprocessing (do not use with Server Datasets)",
                        c("Preprocess", "sc-ATAC-seq")
                    ),

                    conditionalPanel(
                        condition = "input.preprocess_method == 'Preprocess'",
                        ns = ns,
                        actionButton(
                            ns("scanpy"),
                            "See Scanpy's Documentation for more info",
                            class = "longbtn"
                        ),
                        splitLayout(
                            textInput(
                                ns("filter_cells_min"),
                                "Filter cells (min genes)",
                                value = 200
                            ),
                            textInput(
                                ns("filter_cells_max"),
                                "Filter cells (max genes)",
                                value = 3000
                            )
                        ),
                        splitLayout(
                            textInput(
                                ns("filter_genes_min"),
                                "Filter genes (min cells)",
                                value = 3
                            ),
                            textInput(
                                ns("filter_genes_max"),
                                "Filter genes (max cells)",
                                value = -1
                            )
                        ),
                        splitLayout(
                            textInput(
                                ns("normalize_total"),
                                "Normalize total sum",
                                value = 10000
                            ),
                            selectInput(
                                ns("exclude_highly_expressed"),
                                "Exclude highly expressed",
                                c("No", "Yes")
                            )
                        ),
                        div(
                            class = "div_var_genes",
                            tags$p("Highly variable genes (min mean | max mean | min disp)")
                        ),
                        splitLayout(
                            id = "split_layout_var_genes",
                            textInput(
                                ns("high_min_mean"),
                                "",
                                value = 0.0125
                            ),
                            textInput(
                                ns("high_max_mean"),
                                "",
                                value = 3
                            ),
                            textInput(
                                ns("high_min_disp"),
                                "",
                                value = 0.5
                            )
                        ),
                        splitLayout(
                            selectInput(
                                ns("apply_log"),
                                "Logarithmize data",
                                c("Yes", "No")
                            ),
                            textInput(
                                ns("scale"),
                                "Scale max value",
                                value = 10
                            )
                        )
                    ),

                    conditionalPanel(
                        condition = "input.preprocess_method == 'sc-ATAC-seq'",
                        ns = ns,
                        selectInput(
                            ns("atac_operation"),
                            "Operation to apply to bins",
                            c("Sum", "Average")
                        ) %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Whether bin values corresponding ",
                                                "to a gene should be added or averaged."),
                                            placement= "bottom"
                                        )
                                ),
                        splitLayout(
                            textInput(ns("atac_extend"), "Extend", value='5x')
                            %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Can be integer or integer ending in x. ",
                                            "If integer, will expand gene location ",
                                            "interval by that many base pairs. If integer ",
                                            "ending in x, will expand by gene length ",
                                            "times the given integer."),
                                        placement= "bottom"
                                    )
                            ),
                            textInput(ns("atac_max_extend"), "Max Extend", value='50000')
                            %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Same format as Extend.",
                                            "Serves as a threshold for Extend. ",
                                            "This can be useful, for e.g., when ",
                                            "Extend = 5x is specified in terms of ",
                                            "gene length, while Max Extend = 50000 ",
                                            "can take a fixed value."),
                                        placement= "bottom"
                                    )
                            )
                        ),
                        selectInput(
                            ns("atac_strand"),
                            "Consider strand orientation",
                            c("Yes", "No")
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Choose different extension parameters ",
                                            "for the end of the interval that corresponds ",
                                            "to the opposite strand direction of the gene."),
                                        placement= "bottom"
                                    )
                            ),
                        conditionalPanel(
                            condition = "input.atac_strand == 'Yes'",
                            ns = ns,
                            splitLayout(
                                textInput(ns("atac_op_extend"), "Op Extend", value='1x')
                                %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Same format as Extend."),
                                            placement= "bottom"
                                        )
                                ),
                                textInput(ns("atac_max_op_extend"), "Max Op Extend", value='10000')
                                %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Same format as Extend."),
                                            placement= "bottom"
                                        )
                                )
                            )
                        )
                    ),

                    actionButton(
                        ns("runpreprocessbtn"),
                        "Run",
                        class = "sidebtn longbtn")))
        )
    )
)}
