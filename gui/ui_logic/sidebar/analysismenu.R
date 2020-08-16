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
                    
                        selectInput(
                            ns("color"),
                            "View gene expression:",
                            "Clusters"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "View the expression level of the selected gene in the plot. ",
                                        placement= "bottom"
                                    )
                            )
                        ,
                        uiOutput(ns("threshold_slider"))
                        
                        
                ),
            ),

            div(
                class = "div_step div_n_genes",
                list(
                    selectInput(
                        ns("test_method"),
                        "Choose test method",
                        choices = options$de
                    ),
                    sliderInput(
                        ns("mark_markers_n"),
                        "Select number of genes",
                        min = 5, max = 500, value = 50
                    ) %>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Select the number of genes used in the analysis.",
                                    placement= "bottom"
                                )
                        ),
                    splitLayout(
                        cellWidths = c("50%", "50%"),
                        textInput(
                            ns("mark_alpha"),
                            label = "alpha",
                            value = defaults$mark_alpha
                        ),
                        selectInput(
                            ns("mark_correction"),
                            "Correction",
                            choices = options$correction
                        )
                    ),
                    splitLayout(
                        cellWidths = c("50%", "50%"),
                        selectInput(
                            ns("subset1"),
                            "Choose Subset 1",
                            choices = c("None")
                        ) ,
                        # shiny::icon("info-circle") %>%
                        #     bs_embed_tooltip(
                        #         "Subset to run the analysis for against Subset 2.",
                        #         placement= "bottom"
                        #     ),
                        selectInput(
                            ns("subset2"),
                            "Choose Subset 2",
                            choices = c("None")
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
