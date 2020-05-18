library(shiny)

source("gui/options.R")

configmenu <- menuItem(
    "Configuration",
    id = "configbtn",
    icon = icon("cog"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            actionButton(
                "runconfigbtn", "Run with current configuration",
                class="sidebtn longbtn"
            ),

            splitLayout(
                cellWidths = c("60%", "40%"),
                selectInput("dim_method",
                            "Dimensionality reduction:",
                            choices = options$dim),
                textInput(inputId = "dim_n_components",
                          label = "# of components",
                          value = defaults$dim)
            ),

            splitLayout(
                cellWidths = c("60%", "40%"),
                selectInput("clu_method", "Clustering method:",
                            choices = options$clu),
                textInput(inputId = "clu_n_clusters",
                          label = "# of clusters",
                          value = defaults$clu)
            ),

            selectInput(
                "eval_method", "Clustering evaluation method:",
                choices = options$eval
            ),

            splitLayout(
                cellWidths = c("34%", "33%", "33%"),
                textInput(inputId = "mark_alpha",
                          label = "alpha",
                          value = defaults$mark_alpha),
                textInput(inputId = "mark_markers_n",
                          label = "# of markers",
                          value = defaults$mark_no),
                selectInput("mark_correction",
                            "Correction",
                            choices = options$correction)
            ),

            splitLayout(
                cellWidths = c("50%", "50%"),
                selectInput("con_convention",
                            "Converter convention",
                            choices = options$converter),
                textInput(inputId = "con_path",
                          label = "Path (dont use)",
                          value = '')
            ),

            splitLayout(
                cellWidths = c("50%", "50%"),
                selectInput("ide_tissue",
                            "Converter convention",
                            choices = options$tissue),
                textInput(inputId = "ide_path",
                          label = "Path (dont use)",
                          value = "")
            ),

            selectInput(
                "vis_method",
                "Visualization method:",
                choices = options$vis
            ),

            actionButton(
                "ssclurun", "Run constrained clustering",
                class="sidebtn longbtn"
            ),

            splitLayout(
                cellWidths = c("60%", "40%"),
                selectInput(
                    "ssc_method",
                    "Constrained clustering:",
                    choices = options$ssclu
                ),
                textInput(inputId = "savedClusters",
                          label = "Clusters to preserve",
                          value = "")
            )
        )
    )
)
