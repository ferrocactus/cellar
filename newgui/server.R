source("newgui/server_logic/gray_widgets.R")
source("newgui/server_logic/upload_file.R")
source("newgui/server_logic/selected_dataset.R")
source("newgui/server_logic/pipe_cluster.R")
source("newgui/server_logic/gene_card.R")

server <- shinyServer(function(input, output, session) {
    callModule(upload_file, id = "datasetmenu")
    callModule(gray_widgets, id = "configmenu")
    callModule(gene_card, id = "analysismenu")

    selDataset <- callModule(selected_dataset, id = "datasetmenu") # reactive
    pipe <- callModule(cluster_run, id = "configmenu", selDataset = selDataset)
})
