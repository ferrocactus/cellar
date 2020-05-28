library(shiny)
library(shinydashboard)
library(plotly)
library(shinyjs)

source("gui/ui_logic/sidebar/datasetmenu.R") #datasetmenu
source("gui/ui_logic/sidebar/configmenu.R") #configmenu
source("gui/ui_logic/sidebar/selectionmenu.R") #mainmenu
source("gui/ui_logic/sidebar/analysismenu.R") #analysismenu
source("gui/ui_logic/sidebar/downloadmenu.R") #downloadmenu
source("gui/ui_logic/body/plots.R") #plots
source("gui/ui_logic/body/analysis.R") #analysis

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
        datasetmenu(id='ns'),
        configmenu(id='ns'),
        selectionmenu(id='ns'),
        analysismenu(id='ns'),
        downloadmenu(id='ns')
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("gui/ui_logic/styles/style.css")),
    #tags$script(src = "gui/ui_logic/sidebar/anchor.js"),
    plots(id="ns"),
    analysis(id="ns")
)

ui <- dashboardPage(header, sidebar, body, skin='green')
