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

source("newgui/server_logic/go_analysis.R")
source("newgui/server_logic/kegg_analysis.R")

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

    # Functionality
    # We are using the same namespace for everything called "ns".
    # This is not ideal, but a lot of our functions deal with widgets that
    # belong to different ui components, so we are using one namespace for all.
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
    callModule(go_analysis, id = "ns", deGenes = deGenes, pipe = pipe)
    callModule(kegg_analysis, id = "ns", deGenes = deGenes, pipe = pipe)
})
