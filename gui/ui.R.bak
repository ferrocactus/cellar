library(shiny)
library(plotly)

ui <- shinyUI(pageWithSidebar(
  headerPanel("Cellar"),

# Sidebar with a slider input for number of observations
  sidebarPanel(
    shiny::tags$head(shiny::tags$style(shiny::HTML(
      "#log {
          font-size: 15px;
          height: 200px;
          margin-top: 10px;
          overflow: auto;
          padding: 10px;
          background-color: white;
      }
      .select {
          display: inline-block;
      }
      .shiny-split-layout > div {
         overflow: visible;
      }"
    ))),

    splitLayout(
      cellWidths = c("75%", "25%"),
      selectInput("dim_method",
                "Choose a dimensionality reduction method:",
                choices = c("PCA", "UMAP", "TSNE")),
      textInput(inputId = "dim_n_components",
                label = "# of components",
                value = 'knee')
    ),

    splitLayout(
      cellWidths = c("50%", "25%", "25%"),
      selectInput("clu_method",
                "Choose a clustering method:",
                choices = c("KMeans", "KMedoids", "Spectral",
                "Agglomerative", "Birch", "DBSCAN")),
      textInput(inputId = "clu_n_clusters",
                label = "# of clusters",
                value = '(3, 5, 1)'),
      textInput(inputId = "clu_n_jobs",
                label = "# of threads",
                value = -1)
    ),

    selectInput(
      "eval_method",
      "Choose a clustering evaluation method:",
      choices = c("Silhouette", "DaviesBouldin", "CalinskiHarabasz")
    ),

    splitLayout(
      cellWidths = c("25%", "25%", "25%", "25%"),
      textInput(inputId = "mark_alpha",
                label = "alpha",
                value = 0.05),
      textInput(inputId = "mark_markers_n",
                label = "# of markers",
                value = 200),
      selectInput("mark_correction",
                "Multitest Correction",
                choices = c("holm-sidak", "bonferroni", "sidak", "holm", "simes-hochberg",
                "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky")),
      textInput(inputId = "mark_n_jobs",
                label = "# of threads",
                value = 1)
    ),

    splitLayout(
      cellWidths = c("50%", "50%"),
      selectInput("con_convention",
                "Converter convention",
                choices = c("id-to-name", "name-to-id")),
      textInput(inputId = "con_path",
                label = "Path (dont use)",
                value = '')
    ),

    splitLayout(
      cellWidths = c("50%", "50%"),
      selectInput("ide_tissue",
                "Converter convention",
                choices = c("all", "Spleen", "Thyroid", "Kidney", "Liver", "Blood",
                "Placenta", "Eye", "Heart", "Embryo", "Skeletal muscle", "Brain")),
      textInput(inputId = "ide_path",
                label = "Path (dont use)",
                value = '')
    ),

    selectInput(
      "vis_method",
      "Choose a visualization method:",
      choices = c("UMAP", "TSNE")
    ),

    selectInput(
      "ssc_method",
      "Choose a constrained clustering method:",
      choices = c("SeededKMeans")
    ),

    actionButton("update", "Run with current configuration"),
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Plot",
                         h3(textOutput("caption")),
                         plotlyOutput("plot")
                )
    )
  )
))