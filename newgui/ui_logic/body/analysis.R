library(shiny)

analysis <- tabsetPanel(
    id = "switcher",
    # tabPanel(
    #     "No selection",
    #     "No selection"
    # ),
    tabPanel(
        "DE",
        tableOutput("genes"),
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
    #tabPanel(
     #   "Top Expressed Genes",
      #   verbatimTextOutput("topgenes")
    #)
)

