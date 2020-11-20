source("gui/ui_logic/sidebar/options.R")

clusteriimenu <- function(id, label="clusteriimenu") {
ns = NS(id)
menuItem(
    "Constrained Clustering",
    id = "clusteriibtn",
    icon = icon("shapes"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_ssclu",
                list(
                    actionButton(
                        ns("ssclurun"),
                        "Run",
                        class="sidebtn longbtn runconfig"
                    ),
                    selectInput(
                        ns("ssc_method"),
                        "Constrained clustering",
                        choices = options$ssclu
                    )%>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    paste0("Semi-supervised clustering takes into account the ",
                                        "existing labels and tries to refine the clusters."),
                                    placement= "bottom"
                                )
                        ),
                    conditionalPanel(
                        condition = "input.ssc_method != 'SeededKMeans'",
                        ns = ns,
                        splitLayout(
                            textInput(
                                ns("saved_clusters"),
                                "Clusters to preserve",
                                placeholder = "1-5, 8"
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Clusters whose label should not be changed by ",
                                                "the algorithm. Can be integers ",
                                                "or a range specified by a dash, e.g., 1-5."),
                                            placement= "bottom"
                                        )
                                ),
                            textInput(
                                ns("n_ss_clusters"),
                                "Number of clusters"
                            )%>%
                                shinyInput_label_embed(
                                    shiny::icon("info-circle") %>%
                                        bs_embed_tooltip(
                                            paste0("Single integer. Must be greated than the number ",
                                            "of preserved clusters."),
                                            placement="bottom"
                                        )
                                )
                        )
                    )
                )
            ),

            div(
                class = "div_step div_merge",
                list(
                    splitLayout(
                        textInput(
                            ns("clusters_to_merge"),
                            "Clusters to merge",
                            placeholder = "1, 2-5"
                        )%>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        paste0("Integer or range specified by a dash, e.g., 1-5."),
                                        placement= "bottom"
                                    )
                            ),
                        actionButton(
                            ns("merge_clusters"),
                            "Merge clusters",
                            class="sidebtn secondcol"
                        )
                    )
                )
            )
        )
    )
)}