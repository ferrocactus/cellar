# File upload logic
library(reticulate)

source_python("gui/server_logic/read_onto.py")

update_label <- function(input, output, session) {
    dic=get_dic()
    observe({
        updateSelectInput(
            session,
            "tissue",
            "Select tissue",
            choices = sort(names(dic))
        )})
    
    observe({
        
        if (as.character(input$tissue)!=""){
            
            updateSelectInput(
                session,
                "newlabels",
                "Select cell type",
                choices = dic[as.character(input$tissue)]
            )
        }


        })


}

