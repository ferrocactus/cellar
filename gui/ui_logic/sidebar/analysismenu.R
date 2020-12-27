source("gui/ui_logic/sidebar/options.R")

analysismenu <- function(id, label="analysismenu") {
ns = NS(id)
menuItem(
    "Analysis",
    id = "analysisbtn",
    icon = icon("chart-bar"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_color",
                list(
                    splitLayout(
                        selectInput(
                            ns("color"),
                            "View gene expression:",
                            "Clusters"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "View the expression level of the selected gene. ",
                                        placement= "bottom"
                                    )
                            ),
                        selectInput(
                            ns("color2"),
                            "Gene co-expression:",
                            "None"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "View co-expression levels of the two selected genes. ",
                                        placement= "bottom"
                                    )
                            )
                    )
                ),
            ),
            div(
                class = "div_step div_n_genes",
                list(
                    selectInput(
                        ns("test_method"),
                        "Choose test method",
                        choices = options$de
                    ) %>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Select the test to use for finding DE genes.",
                                    placement= "bottom"
                                )
                        ),
                    splitLayout(
                        sliderInput(
                            ns("mark_markers_n"),
                            "Number of genes",
                            min = 5, max = 500, value = 50
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Select the maximum number of DE genes to return. ",
                                                "Note: These genes will all be used to perform enrichment."),
                                        placement= "bottom"
                                    )
                            ),
                        sliderInput(
                            ns("mark_alpha"),
                            "Alpha",
                            min=0.01, max=0.5, step=0.02,
                            value=0.05
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "Select the significance level.",
                                        placement= "bottom"
                                    )
                            )
                    ),
                    # splitLayout(
                    #     cellWidths = c("50%", "50%"),
                    #     textInput(
                    #         ns("mark_alpha"),
                    #         label = "alpha",
                    #         value = defaults$mark_alpha
                    #     ),
                    #     selectInput(
                    #         ns("mark_correction"),
                    #         "Correction",
                    #         choices = options$correction
                    #     )
                    # ),
                    splitLayout(
                        cellWidths = c("50%", "50%"),
                        selectInput(
                            ns("subset1"),
                            "Choose Subset 1",
                            choices = c("All \\ Subset 2")
                        ),
                        # shiny::icon("info-circle") %>%
                        #     bs_embed_tooltip(
                        #         "Subset to run the analysis for against Subset 2.",
                        #         placement= "bottom"
                        #     ),
                        selectInput(
                            ns("subset2"),
                            "Choose Subset 2",
                            choices = c("All \\ Subset 1")
                        )
                        # shiny::icon("info-circle") %>%
                        #     bs_embed_tooltip(
                        #         "If None, will consider all cells not in Subset 1.",
                        #         placement= "bottom"
                        #     )

                    ),
                    actionButton(
                        ns("getdegenes"),
                        "Run Test",
                        class="longbtn",
                        onclick=sprintf("Shiny.onInputChange('%s', 'reset')", ns("select_button")),
                    )
                )
            ),

            div(
                class = "div_step div_genecard",
                list(
                    uiOutput("genecard"),
                    htmlOutput("inc"),
                    splitLayout(
                        textInput(
                            ns("searchgene"),
                            "Search Gene card",
                            placeholder = "Enter gene"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "Gene name to search for.",
                                        placement= "bottom"
                                    )
                            ),
                        actionButton(
                            ns("search"),
                            "Search Card",
                            class = "secondcol"
                        )
                    )
                )
            ),
            div(
                class = "div_step div_genecard",
                list(
                    textInput(
                        ns("n_neighbors"),
                        label = "Neighbors",
                        value = 20
                    )%>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "K nearest neighbors for calculating uncertainty.",
                                    placement= "bottom"
                                )
                        )
                )
            )


        )
    )
)}
