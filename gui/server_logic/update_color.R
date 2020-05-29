# select color value when the plot is updated

update_color <- function(input, output, session, pipe, plotObj,replot) {
  observe({ if (replot() > 0) {
    
    updateSelectInput(session = session,
                      inputId = "color",
                      label = "Select colour value:",
                      choices = c("cluster", as.character(pipe()$col_ids)),
                      selected = NULL)
    
    #replot(0)
  }})
}
