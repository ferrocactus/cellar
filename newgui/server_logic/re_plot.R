re_plot <- function(input, output, session, replot, pipe) {
    observe({
        if (replot() > 0) {
            withProgress(message = "Making plot", value = 0, {
                incProgress(1/4, detail = paste("Step: Rendering plot"))
                output$plot <- renderPlotly({
                    # TODO fix colors and key
                    if (input$color == "cluster") {
                        plotcols = as.factor(pipe()$labels)
                    } else {
                        plotcols = pipe()$labels
                    }

                    plot_ly(
                        x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                        text = ~paste("Label: ", as.factor(pipe()$labels)),
                        color = pipe()$labels,
                        key = plotcols,
                        type = 'scatter',
                        mode = 'markers'
                    ) %>% layout(dragmode = "lasso",
                                title = paste("Value of ", input$color, sep=""))
                })
            })
        }
    })
}
