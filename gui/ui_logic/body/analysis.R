library(shiny)

analysis <- function(id, label="analysis") {
ns = NS(id)
tabsetPanel(
    id = ns("switcher"),
    tabPanel(
        "DE",
        tableOutput(ns("genes")),
        uiOutput(ns("DEbuttons"))
    ),
    tabPanel(
        "Gene Ontology",
        tableOutput(ns("GeneOntology")),
        downloadButton(ns("downloadGO"), "Download GO table")
    ),
    tabPanel(
        "KEGG",
        tableOutput(ns("KEGG")),
        downloadButton(("downloadKEGG"), "Download KEGG table")
    ),
    tabPanel(
        "Markers Intersect",
        tableOutput(ns("Markers")),
        downloadButton(("downloadMKS"), "Download Markers intersect table")
    ),
    tabPanel(
        "MSigDB C2",
        tableOutput(ns("Msigdb")),
        downloadButton(ns("downloadMSIG"), "Download MsigDB enrichment table")
    ),
    tabPanel(
        "User Markers",
         fileInput(ns("markjson"), "Choose markers JSON",multiple = FALSE,accept = c(".json")),
         tableOutput(ns("upmarkers")),
         downloadButton(ns("downloadusermks"), "Download user markers table")
    )
    #tabPanel(
     #   "Top Expressed Genes",
      #   verbatimTextOutput("topgenes")
    #)
)}
