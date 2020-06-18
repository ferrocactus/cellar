source("gui/server_logic/misc.R")
source("gui/server_logic/upload_file.R")
source("gui/server_logic/update_sets.R")
source("gui/server_logic/theme.R")
source("gui/server_logic/singler.R")

source("gui/server_logic/pipe_cluster.R")
source("gui/server_logic/pipe_de.R")
source("gui/server_logic/pipe_align.R")
source("gui/server_logic/plot.R")
source("gui/server_logic/re_mark.R")
source("gui/server_logic/re_label.R")
source("gui/server_logic/analysis_body.R")
source("gui/server_logic/update_label.R")
source("gui/server_logic/save_session.R")
source("gui/server_logic/notifications.R")
source("gui/server_logic/dataset_reset.R")
source("gui/server_logic/download_cells.R")

source("gui/server_logic/dataset.R")
source("gui/server_logic/clustering.R")

# Aliases
is_active <- cellar$utils$r_helpers$is_active
has_key <- cellar$utils$r_helpers$has_key
get_labels <- cellar$utils$r_helpers$get_labels
get_emb_2d <- cellar$utils$r_helpers$get_emb_2d
get_label_names <- cellar$utils$r_helpers$get_label_names
load_file <- cellar$utils$read$load_file
store_subset <- cellar$utils$tools$store_subset
update_subset_label <- cellar$utils$tools$update_subset_label

server <- shinyServer(function(input, output, session) {
    # All variables that need to be used across different modules
    # should be defined here
    adata <- reactiveVal(0)
    selDataset <- reactiveVal("")

    selDatasetAlign <- reactiveVal("")
    pipe <- reactiveVal(0)
    pipeAlign <- reactiveVal(0)
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

    # Dataset menu
    callModule(dataset, id = "ns", adata = adata, selDataset = selDataset,
               selDatasetAlign = selDatasetAlign, fullreset = fullreset)
    callModule(upload_file, id = "ns")

    # Clustering menu
    callModule(cluster, id = "ns", adata = adata, replot = replot,
               reset = reset, relabel = relabel)
    callModule(plot, id = "ns", replot = replot, adata = adata,
               selDataset = selDataset, setNames = setNames, setPts = setPts,
               plotHistory = plotHistory, curPlot = curPlot, reset = reset,
               relabel = relabel)

    # Label Transfer menu
    # callModule(align_run, id = "ns", pipe = pipe, replot = replot,
    #            reset = reset, relabel = relabel,
    #            selDatasetAlign = selDatasetAlign, pipeAlign = pipeAlign)
    # callModule(singler, id = "ns", pipe = pipe, replot = replot,
    #            reset = reset, relabel = relabel,
    #            selDatasetAlign = selDatasetAlign, pipeAlign = pipeAlign)

    # # Selection & Labeling menu
    # callModule(update_sets, id = "ns", setNames = setNames)
    # callModule(update_label, id = "ns", pipe = pipe, labelList = labelList)
    # callModule(re_label, id = "ns", relabel = relabel, pipe = pipe)

    # # Analysis menu
    # callModule(de_run, id = "ns", pipe = pipe, remark = remark,
    #            setNames = setNames, setPts = setPts)
    # callModule(re_mark, id = "ns", remark = remark, pipe = pipe,
    #            deGenes = deGenes, deButtons = deButtons,
    #            rebutton = rebutton, replot = replot)
    # callModule(analysis_body, id = "ns", deGenes = deGenes, pipe = pipe)

    # # Appearance menu
    # callModule(theme, id = "ns", retheme = retheme)

    # # Export/Import menu
    # callModule(save_session, id = "ns", pipe = pipe, setNames = setNames,
    #            setPts = setPts, deGenes = deGenes, selDataset = selDataset,
    #            plotHistory = plotHistory, curPlot = curPlot, replot = replot,
    #            remark = remark, labelList = labelList, relabel = relabel)
    # callModule(download_cells, id = "ns", setNames, setPts,labelList, pipe)

    # # Miscellaneous
    callModule(misc, id = "ns")
    callModule(dataset_reset, id = "ns", reset = reset, setNames = setNames,
               setPts = setPts, labelList = labelList, deButtons = deButtons,
               deGenes = deGenes, adata = adata, fullreset = fullreset,
               curPlot = curPlot, plotHistory = plotHistory)
})
