analysis <- function(id, label="analysis") {
ns = NS(id)
tabsetPanel(
    id = ns("switcher"),
    tabPanel(
        "DE",
        downloadButton(ns("downloadDE"), "Download DE genes"),
        uiOutput(ns("titleDE")),
        DT::dataTableOutput(ns("DEtable")),
    ),
    tabPanel(
      "Heat Map",
      uiOutput(ns("titleheatmap")),
      sliderInput(
        ns("heat_height"),
        "Select heat map height",
        min = 400, max = 800, value = 600
      ),
      plotOutput(ns("heatmap"), height="100%")
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
        uiOutput(ns("titleCellType")),
        downloadButton(ns("downloadCellType"), "Download Cell Type table"),
        DT::dataTableOutput(ns("CellTypetable"))
    ),
    tabPanel(
        "Disease",
        uiOutput(ns("titleDisease")),
        downloadButton(ns("downloadDisease"), "Download Disease table"),
        DT::dataTableOutput(ns("Diseasetable"))
    ),
    tabPanel(
        "User Cell Type",
        downloadButton(ns("downloadUCellType"), "Download Cell Type table"),
        uiOutput(ns("titleUCellType")),
        fileInput(ns("markjson"), "Choose markers JSON",
                  multiple = FALSE, accept = c(".json")),
        DT::dataTableOutput(ns("UCellTypetable"))
    )
)}
