analysis <- function(id, label="analysis") {
ns = NS(id)
tabsetPanel(
    id = ns("switcher"),
    tabPanel(
        "DE",
        conditionalPanel(
            'output.DEtable',
            ns = ns,
            downloadButton(ns("downloadDE"),
                "Export DE genes", class="download_btn_analysis"),
        ),
        uiOutput(ns("titleDE")),
        DT::dataTableOutput(ns("DEtable")),
    ),
    
    tabPanel(
      "Violin Plot",
      sliderInput(
        ns("violin_t"),
        label="Violin plot gene expression thresholds",
        min=0,max=10,
        value=c(0, 10),
        step=0.01
      ),
      uiOutput(ns("titleviolin")),
      plotlyOutput(ns("violin"))
    ),
    
    tabPanel(
      "Heat Map",
      uiOutput(ns("titleheatmap")),
      conditionalPanel(
            'output.heatmap',
            ns = ns,
            sliderInput(
                ns("heat_height"),
                "Select heat map height",
                min = 400, max = 800, value = 600
            )
      ),
      plotOutput(ns("heatmap"), height="100%")
    ),
    tabPanel(
        "Gene Ontology",
        conditionalPanel(
            'output.GOtable',
            ns = ns,
            downloadButton(ns("downloadGO"),
                "Export GO table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleONTO")),
        DT::dataTableOutput(ns("GOtable"))
    ),
    tabPanel(
        "KEGG",
        conditionalPanel(
            'output.KEGGtable',
            ns = ns,
            downloadButton(ns("downloadKEGG"),
                "Export KEGG table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleKEGG")),
        DT::dataTableOutput(ns("KEGGtable"))
    ),
    tabPanel(
        "MSigDB C2",
        conditionalPanel(
            'output.MSIGDBtable',
            ns = ns,
            downloadButton(ns("downloadMSIGDB"),
                "Export MSIGDB table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleMSIG")),
        DT::dataTableOutput(ns("MSIGDBtable"))
    ),
    tabPanel(
        "Disease",
        conditionalPanel(
            'output.Diseasetable',
            ns = ns,
            downloadButton(ns("downloadDisease"),
                "Export Disease table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleDisease")),
        DT::dataTableOutput(ns("Diseasetable"))
    ),
    tabPanel(
        "Cell Type",
        conditionalPanel(
            'output.CellTypetable',
            ns = ns,
            downloadButton(ns("downloadCellType"),
                "Export Cell Type table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleCellType")),
        DT::dataTableOutput(ns("CellTypetable"))
    ),
    tabPanel(
        "External Cell Type",
        conditionalPanel(
            'output.UCellTypetable',
            ns = ns,
            downloadButton(ns("downloadUCellType"),
                "Export Cell Type table", class="download_btn_analysis")
        ),
        uiOutput(ns("titleUCellType")),
        fileInput(ns("markjson"), "Choose markers JSON",
                  multiple = FALSE, accept = c(".json")),
        DT::dataTableOutput(ns("UCellTypetable"))
    )
)}
