source("gui/server_logic/gray_widgets.R")
source("gui/server_logic/upload_file.R")
source("gui/server_logic/selected_dataset.R")
source("gui/server_logic/gene_card.R")
source("gui/server_logic/update_sets.R")
source("gui/server_logic/theme.R")

source("gui/server_logic/load_dataset.R")
source("gui/server_logic/pipe_cluster.R")
source("gui/server_logic/pipe_de.R")
source("gui/server_logic/plot.R")
source("gui/server_logic/re_mark.R")
source("gui/server_logic/re_label.R")
source("gui/server_logic/analysis.R")
source("gui/server_logic/update_label.R")
source("gui/server_logic/analysis_markers.R")
source("gui/server_logic/analysis_disease.R")
source("gui/server_logic/analysis_markers_user.R")
source("gui/server_logic/save_session.R")
source("gui/server_logic/notifications.R")
source("gui/server_logic/dataset_reset.R")
source("gui/server_logic/download_cells.R")

server <- shinyServer(function(input, output, session) {
    # All variables that need to be used across different modules
    # should be defined here
    selDataset <- reactiveVal("")
    pipe <- reactiveVal(0)
    setNames <- reactiveVal(c("None")) # triggers update_sets on change
    setPts <- reactiveVal(c(NA))
    deButtons <- reactiveVal(c())
    deGenes <- reactiveVal(c())
    replot <- reactiveVal(0) # triggers re_plot on change
    remark <- reactiveVal(0) # triggers re_mark and de_buttons on change
    relabel <- reactiveVal(0)
    rebutton <- reactiveVal(0)
    labelList <- reactiveVal(c())
    plotHistory <- reactiveVal(c())
    curPlot <- reactiveVal(0)
    reset <- reactiveVal(0)
    fullreset <- reactiveVal(0)
    retheme <- reactiveVal(0)

    # Functionality
    # We are using the same namespace for everything called "ns".
    # This is not ideal, but a lot of our functions deal with widgets that
    # belong to different ui components, so we are using one namespace for all.
    notificationModule = callModule(notifications, id = 'ns')
    callModule(upload_file, id = "ns")
    callModule(selected_dataset, id = "ns", selDataset = selDataset)
    callModule(load_dataset, id = "ns", pipe = pipe, selDataset = selDataset,
               setNames = setNames, setPts = setPts, fullreset = fullreset)

    callModule(gray_widgets, id = "ns")
    callModule(gene_card, id = "ns")

    callModule(cluster_run, id = "ns", pipe = pipe, selDataset = selDataset,
               replot = replot, reset = reset, relabel = relabel)
    callModule(de_run, id = "ns", pipe = pipe, remark = remark,
               setNames = setNames, setPts = setPts)
    callModule(update_sets, id = "ns", setNames = setNames)
    callModule(update_label, id = "ns", pipe = pipe, labelList = labelList)
    callModule(re_label, id = "ns", relabel = relabel, pipe = pipe)
    callModule(plot, id = "ns", replot = replot, pipe = pipe,
               selDataset = selDataset, setNames = setNames, setPts = setPts,
               plotHistory = plotHistory, curPlot = curPlot, reset = reset,
               relabel = relabel)
    callModule(re_mark, id = "ns", remark = remark, pipe = pipe,
               deGenes = deGenes, deButtons = deButtons,
               rebutton = rebutton, replot = replot)
    callModule(analysis, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(analysis_markers, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(analysis_disease, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(analysis_markers_user, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(dataset_reset, id = "ns", reset = reset, setNames = setNames,
               setPts = setPts, labelList = labelList, deButtons = deButtons,
               deGenes = deGenes, pipe = pipe, fullreset = fullreset,
               curPlot = curPlot, plotHistory = plotHistory)
    callModule(save_session, id = "ns", pipe = pipe, setNames = setNames,
               setPts = setPts, deGenes = deGenes, selDataset = selDataset,
               plotHistory = plotHistory, curPlot = curPlot, replot = replot,
               remark = remark, labelList = labelList, relabel = relabel)
    callModule(theme, id = "ns", retheme = retheme)
    callModule(download_cells, id = "ns", setNames, setPts,labelList, pipe)
})
