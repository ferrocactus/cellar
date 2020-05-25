library(shiny)
library(shinydashboard)
library(plotly)
library(shinyjs)

source("newgui/ui_logic/sidebar/datasetmenu.R") #datasetmenu
source("newgui/ui_logic/sidebar/configmenu.R") #configmenu
source("newgui/ui_logic/sidebar/selectionmenu.R") #mainmenu
source("newgui/ui_logic/sidebar/analysismenu.R") #analysismenu
source("newgui/ui_logic/sidebar/downloadmenu.R") #downloadmenu
source("newgui/ui_logic/body/plots.R") #plots
source("newgui/ui_logic/body/analysis.R") #plots

#load('gui/obj/cell_ontology')

#cell_ontology_names<-paste(cell_ont_full$id,cell_ont_full$name,sep = " ")

header <- dashboardHeader(
    titleWidth = 400,
    title = "Cellar"
)

sidebar <- dashboardSidebar(
    useShinyjs(),
    width = 400,
    sidebarMenu(
        datasetmenu(id='datasetmenu'),
        configmenu(id='configmenu'),
        selectionmenu(id='selectionmenu'),
        analysismenu(id='analysismenu'),
        downloadmenu(id='downloadmenu')
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("newgui/ui_logic/styles/style.css")),
    tags$script(src = "newgui/ui_logic/sidebar/anchor.js"),
    plots,
    analysis
)

ui <- dashboardPage(header, sidebar, body, skin='green')
