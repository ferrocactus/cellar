library(plotly)

scatterplot <- function(input, output, session, x, y, text, color, key) {
    output$plot <- renderPlotly({
        plot_ly(
            x = x, y = y,
            text = text,
            color = color,
            key = key,
            type = 'scatter',
            mode = 'markers'
        ) %>% layout(dragmode = "lasso",
                    title = paste("Value of ", input$color, sep=""))
    })
}
