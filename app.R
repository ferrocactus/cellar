source("gui/ui.R")
source("gui/server.R")

options(shiny.host="0.0.0.0")
options(shiny.port=23123)
options(shiny.maxRequestSize=3000*1024^2)
options(browser = "false")

shinyApp(ui, server)
