library(shiny)
library(plotly)
library(shinyjs)

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

source("gui/mainpanel.R") # mainpanel
source("gui/configpanel.R") # configpanel

ui <- pageWithSidebar(
    headerPanel("Clustering visualization"),

    sidebarPanel(
        id="sidebarpanel",
        includeCSS("gui/style.css"),
        width = 3,

        actionButton("togglemain", "Main Panel", class="panelhead"),
        mainpanel,

        actionButton("toggleconfig", "Configuration", class="panelhead"),
        shinyjs::hidden(configpanel)
    ),

    # Main panel for displaying outputs
    mainPanel(
        tabsetPanel(
            type = "tabs",
            tabPanel(
                "Main Plot",
                h3(textOutput("caption")),
                plotlyOutput("plot"),
                uiOutput("deinfo"),
                verbatimTextOutput("genes"),
                uiOutput("DEbuttons")
            ),
            tabPanel(
                "Updated Plot",
                verbatimTextOutput("brush"),
                plotlyOutput("Plot2")
            ),
            tabPanel(
                "Gene Ontology",
                verbatimTextOutput("GeneOntology")
            ),
            tabPanel(
                "KEGG",
                verbatimTextOutput("KEGG")
            ),
            tabPanel(
                "Markers Intersect",
                verbatimTextOutput("Markers")
            ),
            tabPanel(
                "MSigDB C2",
                verbatimTextOutput("Msigdb")
            )
        ),
        tags$div(id = "placeholder"),
        h3("Clusters and Intersections"),
        tabsetPanel(id = "switcher", tabPanel("No selection", "No selection"))
    )
)
