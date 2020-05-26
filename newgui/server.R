source("newgui/server_logic/gray_widgets.R")
source("newgui/server_logic/upload_file.R")
source("newgui/server_logic/selected_dataset.R")
source("newgui/server_logic/gene_card.R")
source("newgui/server_logic/lasso_store.R")
source("newgui/server_logic/update_sets.R")
source("newgui/server_logic/pipe_cluster.R")
source("newgui/server_logic/pipe_markers.R")
source("newgui/server_logic/re_plot.R")
source("newgui/server_logic/re_mark.R")

server <- shinyServer(function(input, output, session) {
    # Will be reset only when dataset changes
    selDataset <- reactiveVal("")
    setNames <- reactiveVal(c("None")) # triggers update_sets on change
    setPts <- reactiveVal(c(NA))
    pipe <- reactiveVal(0)
    replot <- reactiveVal(0) # triggers re_plot on change
    remark <- reactiveVal(0) # triggers re_mark and de_buttons on change
    deButtons <- reactiveVal(c())
    deGenes <- reactiveVal(c())

    # Functionality
    callModule(upload_file, id = "ns")
    callModule(gray_widgets, id = "ns")
    callModule(gene_card, id = "ns")
    callModule(lasso_store, id = "ns", setNames = setNames, setPts = setPts)

    callModule(selected_dataset, id = "ns", selDataset = selDataset)
    callModule(cluster_run, id = "ns", pipe = pipe,
               selDataset = selDataset, setNames = setNames,
               setPts = setPts, replot = replot, remark = remark,
               deButtons = deButtons)
    callModule(markers_run, id = "ns", pipe = pipe, remark = remark,
               setNames = setNames, setPts = setPts)
    callModule(update_sets, id = "ns", setNames = setNames)
    callModule(re_plot, id = "ns", replot = replot, pipe = pipe)
    callModule(re_mark, id = "ns", remark = remark, pipe = pipe,
               deGenes = deGenes, deButtons = deButtons)
})
