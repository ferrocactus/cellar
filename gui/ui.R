library(shiny)
library(plotly)

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

source("gui/datasetpanel.R") # datasetpanel
source("gui/mainpanel.R") # mainpanel
source("gui/configpanel.R") # configpanel
load('gui/cell_ontology')

#INCLUDE THIS IN THE LABEL DROPDOWN
cell_ontology_names<-paste(cell_ont_full$id,cell_ont_full$name,sep = " ")


ui <- pageWithSidebar(
    headerPanel("Clustering visualization"),

    sidebarPanel(
        id="sidebarpanel",
        includeCSS("gui/style.css"),
        width = 3,

        actionButton("toggledataset", "Dataset", class="panelhead"),
        datasetpanel,

        actionButton("togglemain", "Main Panel", class="panelhead"),
        shinyjs::hidden(mainpanel),

        actionButton("toggleconfig", "Configuration", class="panelhead"),
        shinyjs::hidden(configpanel)
    ),

    # Main panel for displaying outputs
    mainPanel(
        tabsetPanel(
            type = "tabs",
            id = "tabset",
            tabPanel(
                "Main Plot",
                h3(textOutput("caption")),
                plotlyOutput("plot"),
            ),
            tabPanel(
                "Updated Plot",
                verbatimTextOutput("brush"),
                plotlyOutput("Plot2")
            )
            #tabPanel(
             #   "Top Expressed Genes",
              #   verbatimTextOutput("topgenes")
            #),
        ),
        tags$div(id = "placeholder"),
        h3("Clusters and Intersections"),
        tabsetPanel(id = "switcher",
            tabPanel(
                "No selection",
                "No selection"
            ),
            tabPanel(
                "DE",
                verbatimTextOutput("genes"),
                uiOutput("DEbuttons")
            ),
            tabPanel(
                "Gene Ontology",
                tableOutput("GeneOntology"),
                downloadButton("downloadGO", "Download GO table")
            ),
            tabPanel(
                "KEGG",
                tableOutput("KEGG"),
                downloadButton("downloadKEGG", "Download KEGG table")
            ),
            tabPanel(
                "Markers Intersect",
                tableOutput("Markers"),
                downloadButton("downloadMKS", "Download Markers intersect table")
            ),
            tabPanel(
                "MSigDB C2",
                tableOutput("Msigdb"),
                downloadButton("downloadMSIG", "Download MsigDB enrichment table")
            )
        )
    )
)
