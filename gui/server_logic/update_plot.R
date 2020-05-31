update_plot <- function(input, output, session, pipe, setNames,
                        setPts, newLabels) {
    observeEvent(input$labelupd, {
        if (newLabels() == -1) newLabels(pipe()$labels)
        idx=which(setNames()==as.character(input$subset1_upd))
        keys<-setPts()[[idx]]
        #showNotification(as.character(input$newlabels()))
        lbs = newLabels()
        lbs[keys] <- as.character(input$newlabels)
        newLabels(lbs)
        #newLabels()[keys]<-as.character(input$newlabels)
        showNotification("Updating labels")

        #assign("updated_new_labels", newlabs, envir = env)
        if (input$color == 'Clusters') {
          title = "Clusters"
        } else {
          title = input$labelupd
        }
        output$Plot2 <- renderPlotly({
          plot_ly(
            
            x = pipe()$x_emb_2d[, 1], y = pipe()$x_emb_2d[, 2],
            text = ~paste("label: ", as.factor(newLabels())),
            color = as.factor(newLabels()),
            type = 'scatter',
            mode = 'markers'
          ) %>% layout(dragmode = "lasso",
                       title = title)
        })
      })
}
