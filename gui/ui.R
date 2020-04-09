library(shiny)
library(plotly)
library(shinyjs)

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

ui <- pageWithSidebar(
  # App title
  headerPanel("Clustering visualization"),
  ######### Sidebar ###############
  sidebarPanel(
    ############## Style ################
    style="padding: 0;
          border-radius: 0;",

    shiny::tags$head(
      shiny::tags$style(
        shiny::HTML(
          "
          #log {
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
          }
          .sidebtn {
            display: block;
            width: 100%;
          }
          .panelhead {
            border: none;
            background-color: #d8d8d8;
            left: 0;
            border-radius: 0;
            font-size: 16px;
            font-weight: bold;
            height: 50px;
            border-radius: 0;
            width: 100%;
            display: block;
          }
          .panelbody {
            margin-left: 20px;
            margin-right: 20px;
            margin-top: 10px;
            margin-bottom: 10px;
          }
          "
    ))),
    ############# End Style #############
    width = 3,

    actionButton("togglemain", "Main Panel", class="panelhead"),
    #a(id = "togglemain", "Main Panel", href = "#"),

    ########################### Main Block #######################################
    div(id = "mainpanel", class = "panelbody",
      # Input: Selector for the gene to be plotted
      selectInput("color", "Select colour value:", "cluster"),

      sliderInput("nogenes", "Select number of genes",
                    min = 1, max = 500, value = 10),

      splitLayout(
        textInput("text", "New label", placeholder = "Enter label to add..."),
        actionButton("labeladd", "Add", style="margin-top:25px", class="sidebtn")
      ),

      splitLayout(
        selectInput("newlabels", "Select labels", choices = 0),
        actionButton("labelupd", "Update Labels", style="margin-top:25px", class="sidebtn")
      ),

      actionButton("getdegenes", "Get DE genes", class="sidebtn"),

      uiOutput("genecard"),
      htmlOutput("inc"),
      splitLayout(
        textInput("searchgene", "Search Gene card", placeholder = "Enter gene..."),
        actionButton("search", "Search Card", style="margin-top:25px", class="sidebtn")
      ),

      fileInput("file1", "Choose CSV File",
                multiple = FALSE,

                accept = c("text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"))
    ),
    ####################### End of Main Block #####################################


    ######################## Config Block #########################################
    actionButton("toggleconfig", "Configuration", class="panelhead", class="sidebtn"),

    shinyjs::hidden(
      div(id = "configuration", class = "panelbody",
        actionButton("reset", "Run with current configuration",
                    class="sidebtn",
                    style="margin-bottom: 20px;"),

        splitLayout(
          cellWidths = c("60%", "40%"),
          selectInput("dim_method",
                      "Dimensionality reduction method:",
                      choices = c("PCA", "UMAP", "TSNE")),
          textInput(inputId = "dim_n_components",
                    label = "# of components",
                    value = 'knee')
        ),

        splitLayout(
          cellWidths = c("60%", "40%"),
          selectInput("clu_method",
                      "Clustering method:",
                      choices = c("KMeans", "KMedoids", "Spectral",
                              "Agglomerative", "Birch", "DBSCAN",
                              "GaussianMixture")),

          textInput(inputId = "clu_n_clusters",
                    label = "# of clusters",
                    value = '(3, 5, 1)')
        ),

        selectInput(
          "eval_method",
          "Clustering evaluation method:",
          choices = c("Silhouette", "DaviesBouldin", "CalinskiHarabasz")
        ),

        splitLayout(
          cellWidths = c("34%", "33%", "33%"),

          textInput(inputId = "mark_alpha",
                    label = "alpha",
                    value = 0.05),

          textInput(inputId = "mark_markers_n",
                    label = "# of markers",
                    value = 200),

          selectInput("mark_correction",
                      "Multitest Correction",
                      choices = c("holm-sidak", "bonferroni", "sidak", "holm", "simes-hochberg",
                                  "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"))
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

        selectInput(
          "dataset",
          "Choose a dataset:",
          choices = list.files("datasets")
        ),
        useShinyjs(), # Include shinyjs in the UI
        extendShinyjs(text = jsResetCode) # Add the js code to the page
      )
    )
    ###################### End of Config Block ###############################
  ), #end of side bar panel

# Main panel for displaying outputs
  mainPanel(
    tabsetPanel(type = "tabs",

                tabPanel("Main Plot",
                         h3(textOutput("caption")),
                         plotlyOutput("plot"),
                         uiOutput("deinfo"),
                         verbatimTextOutput("genes"),
                         uiOutput("DEbuttons")
                        ),
                tabPanel("Updated Plot",
                         verbatimTextOutput("brush"),
                         plotlyOutput("Plot2")
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
                        )
    ),

    tags$div(id = "placeholder"),
    h3("Clusters and Intersections"),
    tabsetPanel(id = "switcher", tabPanel("No selection", "No selection")
# tabPanel("0", "Cluster0 Genes",tags$div(id = 'placeholder1')),
# tabPanel("1", "Cluster1 Genes",tags$div(id = 'placeholder2')),
# tabPanel("2", "Cluster2 Genes",tags$div(id = 'placeholder3')),
# tabPanel("3", "Cluster3 Genes",tags$div(id = 'placeholder4'))

     ), #end of tabsetpanel
  )
)
