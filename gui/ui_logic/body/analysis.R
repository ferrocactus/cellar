analysis <- function(id, label="analysis") {
ns = NS(id)
tabsetPanel(
    id = ns("switcher"),
    tabPanel(
        "DE",
        downloadButton(ns("downloadDE"), "Download DE genes"),
        uiOutput(ns("titleDE")),
        DT::dataTableOutput(ns("DEtable")),
        #uiOutput(ns("DEbuttons"))
    ),
    tabPanel(
      "Heat Map",
      uiOutput(ns("titleheatmap")),
      plotOutput(ns("heatmap"), height="600px")
    ),
    tabPanel(
        "Gene Ontology",
        downloadButton(ns("downloadGO"), "Download GO table"),
        uiOutput(ns("titleONTO")),
        DT::dataTableOutput(ns("GOtable"))
    ),
    tabPanel(
        "KEGG",
        downloadButton(ns("downloadKEGG"), "Download KEGG table"),
        uiOutput(ns("titleKEGG")),
        DT::dataTableOutput(ns("KEGGtable"))
    ),
    tabPanel(
        "MSigDB C2",
        downloadButton(ns("downloadMSIGDB"), "Download MSIGDB table"),
        uiOutput(ns("titleMSIG")),
        DT::dataTableOutput(ns("MSIGDBtable"))
    ),
    tabPanel(
        "Cell Type",
        uiOutput(ns("titleM")),
        downloadButton(ns("downloadMKS"), "Download Markers table"),
        DT::dataTableOutput(ns("Markerstable"))
    ),
    tabPanel(
        "Disease",
        uiOutput(ns("titleD")),
        downloadButton(ns("downloadDis"), "Download Markers table"),
        DT::dataTableOutput(ns("disease"))
    ),
    tabPanel(
        "User Markers",
        downloadButton(ns("downloadMKSuser"), "Download Markers table"),
        uiOutput(ns("titleUM")),
        fileInput(ns("markjson"), "Choose markers JSON",
                  multiple = FALSE, accept = c(".json")),
        DT::dataTableOutput(ns("Markerstableuser"))
    )
    #tabPanel(
     #   "Top Expressed Genes",
      #   verbatimTextOutput("topgenes")
    #)
)}
