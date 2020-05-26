lasso_store <- function(input, output, session, setNames, setPts) {
    observeEvent(input$store_lasso, {
        if (as.character(input$newsubset) %in% setNames()){
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
          setPts(c(setPts(),list(keys)))
          setNames(c(setNames(),as.character(input$newsubset)))

        }
        else{
          setPts(c(setPts(),list(keys)))
          setNames(c(setNames(),paste("Newset",as.character(length(setNames())),sep="")))
        }

        showNotification(paste(as.character(cell_count)," cells stored",sep=""))
    })
}
