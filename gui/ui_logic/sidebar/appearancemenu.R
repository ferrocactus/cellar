library(shiny)

appearancemenu <- function(id, label='appearancemenu') {
  ns = NS(id)
  menuItem(
    "Appearance",
    id = "appearancebtn",
    icon = icon("archive"),
    startExpanded = TRUE,
    menuSubItem(
      icon = NULL,
      list(
        
        sliderInput(
          ns("dot_size"),
          "Select dot size",
          min = 1, max = 30, value = 3
        )
        
      )
    )
  )}
