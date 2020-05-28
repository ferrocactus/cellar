source("gui/server_logic/gray_widgets.R")
source("gui/server_logic/upload_file.R")
source("gui/server_logic/selected_dataset.R")
source("gui/server_logic/gene_card.R")
source("gui/server_logic/lasso_store.R")
source("gui/server_logic/update_sets.R")
source("gui/server_logic/add_label.R")

source("gui/server_logic/pipe_cluster.R")
source("gui/server_logic/pipe_de.R")
source("gui/server_logic/re_plot.R")
source("gui/server_logic/re_mark.R")
source("gui/server_logic/analysis.R")
source("gui/server_logic/update_label.R")
source("gui/server_logic/update_plot.R")
source("gui/server_logic/analysis_markers.R")
source("gui/server_logic/save_session.R")

server <- shinyServer(function(input, output, session) {
    # All variables that need to be used across different modules
    # should be defined here
    selDataset <- reactiveVal("")
    setNames <- reactiveVal(c("None")) # triggers update_sets on change
    setPts <- reactiveVal(c(NA))
    pipe <- reactiveVal(0)
    replot <- reactiveVal(0) # triggers re_plot on change
    remark <- reactiveVal(0) # triggers re_mark and de_buttons on change
    deButtons <- reactiveVal(c())
    deGenes <- reactiveVal(c())
    labelList <- reactiveVal(c())
    plotObj <- reactiveVal(0)
    newLabels <- reactiveVal(-1)

    # Functionality
    # We are using the same namespace for everything called "ns".
    # This is not ideal, but a lot of our functions deal with widgets that
    # belong to different ui components, so we are using one namespace for all.
    callModule(upload_file, id = "ns")
    callModule(gray_widgets, id = "ns")
    callModule(gene_card, id = "ns")
    callModule(lasso_store, id = "ns", setNames = setNames, setPts = setPts)
    callModule(add_label, id = "ns", labelList = labelList)

    callModule(selected_dataset, id = "ns", selDataset = selDataset)
    callModule(cluster_run, id = "ns", pipe = pipe,
               selDataset = selDataset, setNames = setNames,
               setPts = setPts, replot = replot, remark = remark,
               deButtons = deButtons, deGenes = deGenes, labelList = labelList)
    callModule(de_run, id = "ns", pipe = pipe, remark = remark,
               setNames = setNames, setPts = setPts)
    callModule(update_sets, id = "ns", setNames = setNames)
    callModule(update_label, id = "ns", pipe = pipe, labelList = labelList)
    callModule(update_plot, id = "ns", pipe = pipe, setNames = setNames,
               setPts = setPts, newLabels = newLabels)
    callModule(re_plot, id = "ns", replot = replot, pipe = pipe,
               plotObj = plotObj, selDataset = selDataset)
    callModule(re_mark, id = "ns", remark = remark, pipe = pipe,
               deGenes = deGenes, deButtons = deButtons)
    callModule(analysis, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(analysis_markers, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(save_session, id = "ns", pipe = pipe, setNames = setNames,
               setPts = setPts, deGenes = deGenes, selDataset = selDataset)
})
