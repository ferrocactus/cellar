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

      uiOutput('violin_t'),
      sliderInput(
        ns("violin_t"),
        label="Violin plot gene expression thresholds",
        min=-1,max=10,
        value=c(4.99, 5.11),
        step=0.01
      ),
      uiOutput(ns("titleviolin")),
      imageOutput(ns("violin")),#,width = '800px', height = "400px",),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
      br(),
       plotlyOutput(ns("zeros"))
    ),

    tabPanel(
        "Heat Map",
        uiOutput(ns("titleheatmap")),
        plotOutput(ns("heatmap"), height="100%"),
        multiInput(
            ns("heatmap_genes"),
            "Select Genes to use in Heatmap",
            choices=c("None")
        ),
        splitLayout(
            sliderInput(
                ns("heat_height"),
                "Heatmap Height",
                min = 400, max = 800, value = 600
            ),
            actionButton(
                ns("build_heatmap"),
                "Generate Heatmap w Selected Genes",
                class="heatmap_btn"
            ),
            actionButton(
                ns("append_de"),
                "Append Current DE Genes",
                class="heatmap_btn"
            ),
            actionButton(
                ns("clear_selected_genes"),
                "Clear Selected Genes",
                class="heatmap_btn"
            )
        )
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
        "User Defined",
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
    ),
    tabPanel(
        "CODEX",
        splitLayout(
            fileInput(ns("codex_upload"), "Upload CODEX data",
                    multiple = FALSE, accept = c(
                        ".tar.gz", "application/tar+gzip", "application/gzip")),
            actionButton(
                ns("codex_generate"),
                "(Re)Generate CODEX Tile"
            )
        ),
        imageOutput(ns('codex_tile')),
        style='height: 1700px'
    )
)}
