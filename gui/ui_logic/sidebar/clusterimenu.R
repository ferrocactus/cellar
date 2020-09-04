source("gui/ui_logic/sidebar/options.R")

clusterimenu <- function(id, label="clusterimenu") {
ns = NS(id)
menuItem(
    "Clustering",
    id = "clusteribtn",
    icon = icon("shapes"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_dim",
                list(
                    actionButton(
                        ns("runconfigbtn"),
                        "Run",
                        class="sidebtn longbtn runconfig"
                    ),
                    selectInput(
                        ns("dim_method"),
                        "Dimensionality reduction",
                        choices = options$dim,
                        selected='Truncated SVD'
                    )%>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Choose a dimension reduction method",
                                    placement= "bottom"
                                )
                        ),
                    conditionalPanel(
                        condition = "input.dim_method != 'Monocle3'",
                        ns = ns,
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
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("If automatic, will try and detect the optimal ",
                                                "number of components using the ankle heuristic of finding the ",
                                                "ankle in the explained variance graph."),
                                            placement= "bottom"
                                        )
                                ),
                            div(
                                class = 'div_dim_n_components',
                                textInput(
                                    ns("dim_n_components"),
                                    label = NULL,
                                    value = defaults$dim
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            "Integer value.",
                                            placement= "bottom"
                                        )
                                )

                            )
                        )
                    )
                )
            ),

            div(
                class = "div_step div_clu",
                list(
                    selectInput(
                        ns("clu_method"),
                        "Clustering",
                        choices = options$clu
                    )%>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Choose a clustering method",
                                    placement= "bottom"
                                )
                        ),
                    conditionalPanel(
                        condition = "input.clu_method == 'Leiden'",
                        ns = ns,
                        splitLayout(
                            textInput(
                                ns("leiden_resolution"),
                                "Resolution",
                                value = 1
                            ),
                            textInput(
                                ns("leiden_neighbors"),
                                "Number of neighbors",
                                value = 15
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = "input.clu_method != 'Leiden'",
                        ns = ns,
                        splitLayout(
                            cellWidths = c("50%", "50%"),
                            textInput(
                                ns("clu_n_clusters"),
                                "Number of clusters",
                                value = defaults$clu
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Can be a single integer, ",
                                                "a list of comma separated values, ",
                                                "or a tuple specifying a range, e.g., (3, 9, 2) will start ",
                                                "at 3 and finish at 9 on increments of 2."),
                                            placement= "bottom"
                                        )
                                ),
                            selectInput(
                                ns("eval_method"),
                                "Evaluation:",
                                choices = options$eval
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Used to determine the best number of clusters ",
                                                "if number of clusters is a list."),
                                            placement= "bottom"
                                        )
                                )
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
                )
            ),

            div(
                class = "div_step div_vis",
                list(
                    selectInput(
                        ns("vis_method"),
                        "Visualization method:",
                        choices = options$vis
                    )
                )
            )
        )
    )
)}
