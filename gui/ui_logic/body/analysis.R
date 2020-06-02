analysis <- function(id, label="analysis") {
ns = NS(id)
tabsetPanel(
    id = ns("switcher"),
    tabPanel(
        "DE",
        downloadButton(ns("downloadDE"), "Download DE genes"),
        tableOutput(ns("DEtable")),
        uiOutput(ns("DEbuttons"))
    ),
    tabPanel(
        "Gene Ontology",
        downloadButton(ns("downloadGO"), "Download GO table"),
        tableOutput(ns("GOtable"))
    ),
    tabPanel(
        "KEGG",
        downloadButton(ns("downloadKEGG"), "Download KEGG table"),
        tableOutput(ns("KEGGtable"))
    ),
    tabPanel(
        "MSigDB C2",
        downloadButton(ns("downloadMSIGDB"), "Download MSIGDB table"),
        tableOutput(ns("MSIGDBtable"))
    ),
    tabPanel(
        "Markers",
        downloadButton(ns("downloadMKS"), "Download Markers table"),
        tableOutput(ns("Markerstable"))
    ),
    tabPanel(
        "User Markers",
        downloadButton(ns("downloadMKSuser"), "Download Markers table"),
        fileInput(ns("markjson"), "Choose markers JSON",
                  multiple = FALSE, accept = c(".json")),
        tableOutput(ns("Markerstableuser"))
    )
    #tabPanel(
     #   "Top Expressed Genes",
      #   verbatimTextOutput("topgenes")
    #)
)}
