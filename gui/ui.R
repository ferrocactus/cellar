library(shiny)
library(plotly)
library(shinyjs)

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

ui <- pageWithSidebar(
  # App title
  headerPanel("Clustering visualization"),

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

                      # Input: Selector for the gene to be plotted
                       selectInput("color", "Select colour value:",
                       "cluster"),

                      sliderInput("nogenes", "Select number of genes",
                                   min = 1, max = 100, value = 10
                      ),

                      selectInput("newlabels","Select labels", choices=0),

                      textInput("text", "New label", value = "Enter label..."),

                      actionButton("labeladd", "Add"),

                      actionButton("labelupd", "Update Labels"),

                      actionButton("getdegenes", "Get DE genes"),

                      width = 2,
                      uiOutput("genecard"),
                      htmlOutput("inc"),
                      textInput("searchgene", "Search Gene card", value = "Enter gene..."),
                      actionButton("search","Search Card"),
                      fileInput("file1", "Choose CSV File",
                               multiple = FALSE,

                               accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
  ), #end of side bar panel

  # Main panel for displaying outputs
  mainPanel(
    tabsetPanel(type = "tabs",

                tabPanel("Main Plot",
                         h3(textOutput("caption")),
                         plotlyOutput("plot")
                        ),
                tabPanel("Updated Plot",
                         verbatimTextOutput("brush"),
                         plotlyOutput("Plot2")
                         ),
                tabPanel("DE genes",
                         #verbatimTextOutput("click"),
                         verbatimTextOutput("genes")
                         ),
                tabPanel("Gene Ontology",
                         verbatimTextOutput("GeneOntology")
                         ),
                tabPanel("KEGG",
                         verbatimTextOutput("KEGG")),
                tabPanel("Markers Intersect",
                         verbatimTextOutput("Markers")
                         ),
                tabPanel("MSigDB C2",
                         verbatimTextOutput("Msigdb")
                        ),
                tabPanel(
                  "Configurations",
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
                              value = 1)
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
                  useShinyjs(),                                           # Include shinyjs in the UI
                  extendShinyjs(text = jsResetCode),                      # Add the js code to the page
                  actionButton("reset", "Run with current configuration")
                  #actionButton("update", "Run with current configuration"),
                )


    ),
     uiOutput("DEbuttons"),
     tags$div(id="placeholder"),
     tabsetPanel(
       id = "switcher",
       tabPanel("No selection", "No selection")
       # tabPanel("0", "Cluster0 Genes",tags$div(id = 'placeholder1')),
       # tabPanel("1", "Cluster1 Genes",tags$div(id = 'placeholder2')),
       # tabPanel("2", "Cluster2 Genes",tags$div(id = 'placeholder3')),
       # tabPanel("3", "Cluster3 Genes",tags$div(id = 'placeholder4'))

     ), #end of tabsetpanel
  )
)


