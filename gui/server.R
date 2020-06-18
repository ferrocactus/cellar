source("gui/server_logic/misc.R")
source("gui/server_logic/upload_file.R")
source("gui/server_logic/theme.R")
source("gui/server_logic/singler.R")

source("gui/server_logic/pipe_cluster.R")
source("gui/server_logic/pipe_de.R")
source("gui/server_logic/pipe_align.R")
source("gui/server_logic/re_mark.R")
source("gui/server_logic/analysis_body.R")
source("gui/server_logic/save_session.R")
source("gui/server_logic/notifications.R")
source("gui/server_logic/dataset_reset.R")
source("gui/server_logic/download_cells.R")

source("gui/server_logic/dataset.R")
source("gui/server_logic/clustering.R")
source("gui/server_logic/plot.R")

source("gui/server_logic/selectionLabeling.R")

# Aliases
is_active <- cellar$utils$r_helpers$is_active
has_key <- cellar$utils$r_helpers$has_key
get_labels <- cellar$utils$r_helpers$get_labels
get_emb_2d <- cellar$utils$r_helpers$get_emb_2d
get_label_names <- cellar$utils$r_helpers$get_label_names
load_file <- cellar$utils$read$load_file
store_subset <- cellar$utils$tools$store_subset
update_subset_label <- cellar$utils$tools$update_subset_label
get_cluster_label_list <- cellar$utils$r_helpers$get_cluster_label_list
get_cluster_name_list <- cellar$utils$r_helpers$get_cluster_name_list
get_unique_labels <- cellar$utils$r_helpers$get_unique_labels
get_subsets <- cellar$utils$r_helpers$get_subsets

server <- shinyServer(function(input, output, session) {
    # All variables that need to be used across different modules
    # should be defined here
    adata <- reactiveVal(0)
    selDataset <- reactiveVal("")
    labelList <- reactiveVal(c())

    selDatasetAlign <- reactiveVal("")
    pipe <- reactiveVal(0)
    pipeAlign <- reactiveVal(0)
    setNames <- reactiveVal(c("None")) # triggers update_sets on change
    deGenes <- reactiveVal(c())
    replot <- reactiveVal(0) # triggers re_plot on change
    remark <- reactiveVal(0) # triggers re_mark and de_buttons on change
    relabel <- reactiveVal(0)
    rebutton <- reactiveVal(0)
    plotHistory <- reactiveVal(c())
    curPlot <- reactiveVal(0)
    reset <- reactiveVal(0)
    fullreset <- reactiveVal(0)
    retheme <- reactiveVal(0)
    resubset <- reactiveVal(0)

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
               reset = reset, relabel = relabel, resubset = resubset)
    callModule(plot, id = "ns", replot = replot, adata = adata,
               selDataset = selDataset, setNames = setNames, setPts = setPts,
               plotHistory = plotHistory, curPlot = curPlot, reset = reset,
               resubset = resubset)

    # Label Transfer menu
    # callModule(align_run, id = "ns", pipe = pipe, replot = replot,
    #            reset = reset, relabel = relabel,
    #            selDatasetAlign = selDatasetAlign, pipeAlign = pipeAlign)
    # callModule(singler, id = "ns", pipe = pipe, replot = replot,
    #            reset = reset, relabel = relabel,
    #            selDatasetAlign = selDatasetAlign, pipeAlign = pipeAlign)

    # # Selection & Labeling menu
    callModule(selectionLabeling, id = "ns", adata = adata,
               labelList = labelList, setNames = setNames,
               resubset = resubset, reset = reset, replot = replot,
               relabel = relabel)

    # # Analysis menu
    # callModule(de_run, id = "ns", pipe = pipe, remark = remark,
    #            setNames = setNames, setPts = setPts)
    # callModule(re_mark, id = "ns", remark = remark, pipe = pipe,
    #            deGenes = deGenes,
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
               labelList = labelList,
               deGenes = deGenes, adata = adata, fullreset = fullreset,
               curPlot = curPlot, plotHistory = plotHistory,
               resubset = resubset)
})
