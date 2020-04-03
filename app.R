source("gui/ui.R")
source("gui/server.R")

options(shiny.host="0.0.0.0")
options(shiny.port=23123)

shinyApp(ui, server)
