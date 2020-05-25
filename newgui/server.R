source("newgui/server_logic/gray_widgets.R")
source("newgui/server_logic/upload_file.R")
source("newgui/server_logic/selected_dataset.R")
source("newgui/server_logic/pipe_cluster.R")
source("newgui/server_logic/gene_card.R")
source("newgui/server_logic/lasso_store.R")

server <- shinyServer(function(input, output, session) {
    callModule(upload_file, id = "ns")
    callModule(gray_widgets, id = "ns")
    callModule(gene_card, id = "ns")
    #callModule(lasso_store, id = "selectionmmenu analysismenu plots")

    selDataset <- callModule(selected_dataset, id = "ns") # reactive
    pipe <- callModule(cluster_run, id = "ns", selDataset = selDataset)
})
