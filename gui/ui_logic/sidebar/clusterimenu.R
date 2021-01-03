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
                    selectInput(
                        ns("dim_method"),
                        "Dimensionality reduction",
                        choices = options$dim,
                        selected='PCA'
                    )%>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    paste0("Choose a dimension reduction method. These ",
                                            "embeddings will be used in the clustering ",
                                            "and visualization steps."),
                                    placement= "bottom"
                                )
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
                            selected = 'pca_manual',
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
                            sliderInput(
                                ns("dim_n_components"),
                                label = NULL,
                                min=2, max=100, step=1,
                                value=40
                            ) %>%
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
            ),

            div(
                class = "div_step div_vis",
                list(
                    selectInput(
                        ns("vis_method"),
                        "Visualization method:",
                        choices = options$vis
                    ) %>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    paste0("Choose a reduction method to further reduce ",
                                            "the data to 2 dimensions for visualization."),
                                    placement= "bottom"
                                )
                        )
                )
            ),

            actionButton(
                ns("reduce_dim_and_vis"),
                "Reduce Dimensions and Visualize",
                class="sidebtn longbtn runconfig"
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
                                    paste0("Choose a clustering method. The ",
                                            "clustering is applied to the reduced ",
                                            "data obtained via the first dimensionality ",
                                            "reduction method selected above (i.e., not ",
                                            "the 2D embeddings)."),
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
                            ) %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("A floating point value. ",
                                            "A higher value results in more clusters. ",
                                            "Consider using 0.1-0.5 for few clusters."),
                                            placement= "bottom"
                                        )
                                ),
                            sliderInput(
                                ns("leiden_neighbors"),
                                "Number of neighbors",
                                min=5, max=100, step=5,
                                value = 15
                            ) %>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Leiden takes a connectivity graph as input. ",
                                                    "We construct this as a neighbors graph. ",
                                                    "Here you can specify the number of neighbors ",
                                                    "to use in the graph."),
                                            placement= "bottom"
                                        )
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
                ),

                actionButton(
                    ns("run_clustering"),
                    "Cluster",
                    class="sidebtn longbtn runconfig"
                ),
            )
        )
    )
)}
