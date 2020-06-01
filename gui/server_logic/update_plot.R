update_plot <- function(input, output, session, pipe, setNames,
                        setPts, newLabels, plotObj) {
    observeEvent(input$labelupd, {
        if (is.null(newLabels())) newLabels(pipe()$labels)
        idx=which(setNames()==as.character(input$subset1_upd))
        keys<-setPts()[[idx]]
        lbs = newLabels()
        lbs[keys] <- as.character(input$newlabels)
        newLabels(lbs)
        showNotification("Updating labels")

        updateCheckboxGroupInput(
            session,
            'saved_clusters',
            label = "Select clusters to preserve",
            choices = sort(unique(newLabels()))
        )

        if (input$color == 'Clusters') {
          title = "Clusters"
        } else {
          title = input$labelupd
        }
        ns <- session$ns
        output$plot <- renderPlotly({
          plotObj(plot_ly(
            x = pipe()$x_emb_2d[, 1], y = pipe()$x_emb_2d[, 2],
            text = ~paste("label: ", as.factor(newLabels())),
            color = as.factor(newLabels()),
            type = 'scatter',
            mode = 'markers',
            source = ns("M")
          ) %>% layout(dragmode = "lasso",
                       title = title))
          return(plotObj())
        })
      })
}
