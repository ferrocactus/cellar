# File upload logic
library(reticulate)

source_python("read_onto.py")

update_label <- function(input, output, session) {
    observe({
        updateSelectInput(
            session,
            "tissue",
            "Select tissue",
            choices = sort(names(dic))
        )})
    dic['']=''
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

