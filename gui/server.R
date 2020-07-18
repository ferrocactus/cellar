source("gui/server_logic/misc.R")
source("gui/server_logic/upload_file.R")
source("gui/server_logic/theme.R")

source("gui/server_logic/analysis_body.R")
source("gui/server_logic/save_session.R")
source("gui/server_logic/notifications.R")
source("gui/server_logic/dataset_reset.R")
source("gui/server_logic/download_cells.R")

source("gui/server_logic/dataset.R")
source("gui/server_logic/clustering.R")
source("gui/server_logic/plot.R")
source("gui/server_logic/align.R")
source("gui/server_logic/selectionLabeling.R")
source("gui/server_logic/de.R")
source("gui/server_logic/uncertainty.R")
cellar <- import("src", convert=FALSE)
source("gui/server_logic/aliases.R")

server <- shinyServer(function(input, output, session) {
    # All variables that need to be used across different modules
    # should be defined here
    adata <- reactiveVal(0)

    selDataset <- reactiveVal("")
    selDatasetAlign <- reactiveVal("")
    labelList <- reactiveVal(c())
    setNames <- reactiveVal(c("None"))
    deGenes <- reactiveVal(0)

    reset <- reactiveVal(0)
    fullreset <- reactiveVal(0)
    replot <- reactiveVal(0)
    remark <- reactiveVal(0)
    relabel <- reactiveVal(0)
    reinfo <- reactiveVal(0)
    resubset <- reactiveVal(0)
    retheme <- reactiveVal(0)

    curPlot <- reactiveVal(0)
    plotHistory <- reactiveVal(c())

    cellNamesTb <- reactiveVal(c())
    infoTb <- reactiveVal(c())

    # We are using the same namespace for everything called "ns".
    notificationModule = callModule(notifications, id = 'ns')

    # Dataset menu
    callModule(dataset, id = "ns", adata = adata, selDataset = selDataset,
               selDatasetAlign = selDatasetAlign, fullreset = fullreset)
    callModule(upload_file, id = "ns")

    # Clustering menu
    callModule(cluster, id = "ns", adata = adata, replot = replot,
               reset = reset, relabel = relabel, resubset = resubset,
               reinfo = reinfo,
               cellNamesTb = cellNamesTb, infoTb = infoTb)
    callModule(plot, id = "ns", replot = replot, adata = adata,
               selDataset = selDataset, setNames = setNames,
               plotHistory = plotHistory, curPlot = curPlot,
               reset = reset, resubset = resubset,
               cellNamesTb = cellNamesTb, infoTb = infoTb,
               reinfo = reinfo, relabel = relabel)

    # Label Transfer menu
    callModule(align, id = "ns", adata = adata,
               selDatasetAlign = selDatasetAlign, replot = replot,
               reset = reset, relabel = relabel, resubset = resubset,
               reinfo = reinfo)

    # Selection & Labeling menu
    callModule(selectionLabeling, id = "ns", adata = adata,
               labelList = labelList, setNames = setNames,
               resubset = resubset, reset = reset, replot = replot,
               relabel = relabel, reinfo = reinfo)

    # Analysis menu
    callModule(differential_e, id = "ns", adata = adata, remark = remark,
               deGenes = deGenes)
    callModule(analysis_body, id = "ns", adata = adata, deGenes = deGenes)
    callModule(uncertain, id = "ns", adata = adata)
    
    # Appearance menu
    callModule(theme, id = "ns", retheme = retheme)

    # Export/Import menu
    callModule(save_session, id = "ns", adata = adata, replot = replot,
               remark = remark, labelList = labelList, relabel = relabel,
               resubset = resubset, fullreset = fullreset, reinfo = reinfo)
    callModule(download_cells, id = "ns", adata = adata)

    # # Miscellaneous
    callModule(misc, id = "ns")
    callModule(dataset_reset, id = "ns", reset = reset, setNames = setNames,
               labelList = labelList, deGenes = deGenes, adata = adata,
               fullreset = fullreset, curPlot = curPlot,
               plotHistory = plotHistory, resubset = resubset)
})
