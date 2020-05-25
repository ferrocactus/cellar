lasso_store <- function(input, output, session) {
    set_name <- c("None")
    sets <- c(NA)

      for (i in 1:length(markers)){
        assign("set_name",c(set_name,(paste("cluster_",as.character(i-1),sep=""))),envir=env)
        keys<-which(expr_data$cluster==(i-1))
        assign("sets",c(sets,list(keys)),envir=env)
      }

    observeEvent(input$store_lasso, {
        if (as.character(input$newsubset) %in% set_name){
          showNotification("Name already exists")
          return()
        }

        if (substr(as.character(input$newsubset),1,7)=="cluster"){
          showNotification("Reserved name, please choose another.")
          return()
        }

        d <- event_data("plotly_selected")
        keys=as.numeric(d$key)
        cell_count=length(d$key)
        #showNotification(as.character(input$newsubset),duration=NULL)
        if (identical(d$key,NULL)==TRUE){
          return()
        }
        if (as.character(input$newsubset)!=""){
          sets<<-c(sets,list(keys))
          set_name<<-c(set_name,as.character(input$newsubset))

        }
        else{
          sets<<-c(sets,list(keys))
          set_name<<-c(set_name,paste("Newset",as.character(length(set_name)),sep=""))
        }

        showNotification(paste(as.character(cell_count)," cells stored",sep=""))

        # select subset1
        updateSelectInput(session = session,
                          inputId = "subset1",
                          label = "Choose Subset 1",
                          choices = set_name,
                          selected = NULL)
        # select subset2
        updateSelectInput(session = session,
                          inputId = "subset2",
                          label = "Choose Subset 2",
                          choices = set_name,
                          selected = NULL)

        updateSelectInput(session = session,
                          inputId = "cell_subset_download",
                          label = "Download labels",
                          choices = c("All", set_name),
                          selected = 'All')

    })
}
