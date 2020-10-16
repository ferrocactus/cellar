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
                        c("Defaults", "sc-ATAC-seq", "Manual")
                    ),
                    conditionalPanel(
                        condition = "input.preprocess_method == 'Manual'",
                        ns = ns,
                        tags$p("Coming soon...")
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
