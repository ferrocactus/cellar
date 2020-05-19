library(shiny)
library(shinydashboard)
library(plotly)
library(shinyjs)

source("gui/datasetmenu.R") #datasetmenu
source("gui/mainmenu.R") #mainmenu
source("gui/configmenu.R") #configmenu

load('gui/cell_ontology')

cell_ontology_names<-paste(cell_ont_full$id,cell_ont_full$name,sep = " ")

header <- dashboardHeader(
    titleWidth = 400,
    title = "Cellar"
)

sidebar <- dashboardSidebar(
    width = 400,
    sidebarMenu(
        datasetmenu,
        mainmenu,
        configmenu
    )
)

body <- dashboardBody(
    tags$head(includeCSS("./gui/custom.css")),
    # Don't scroll to top after clicking anything
    tags$script(HTML(
        "
        $(document).on('click', function(e) {
           e.preventDefault();
        });
        "
    )),
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

ui <- dashboardPage(header, sidebar, body, skin='green')
