# File upload logic
library(reticulate)

source_python("gui/server_logic/read_onto.py")

update_label <- function(input, output, session,pipe,labelList) {
    dic=get_dic()
    observe({
        updateSelectInput(
            session,
            "tissue",
            "Select tissue",
            choices = c(sort(names(dic)),'Clusters','User defined')
        )})
    
    observe({
        print(input$tissue)
        if (as.character(input$tissue)!=""){
            if (as.character(input$tissue)=='Clusters')
            {
                n_clusters=pipe()['n_clusters']
                updateSelectInput(
                    session,
                    "newlabels",
                    "Select cell type",
                    choices = n_clusters
                )
            }
            else if (as.character(input$tissue)=='User defined')
            {
                n_clusters=pipe()['n_clusters']
                updateSelectInput(
                    session,
                    "newlabels",
                    "Select cell type",
                    choices = n_clusters
                )
            }
            else{
                updateSelectInput(
                    session,
                    "newlabels",
                    "Select cell type",
                    choices = dic[as.character(input$tissue)]
                )    
            }
            
        }
        
        
    })


}

