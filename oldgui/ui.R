library(shiny)
library(shinydashboard)
library(plotly)
library(shinyjs)

source("gui/sidebar/datasetmenu.R") #datasetmenu
source("gui/sidebar/configmenu.R") #configmenu
source("gui/sidebar/selectionmenu.R") #mainmenu
source("gui/sidebar/analysismenu.R") #analysismenu
source("gui/sidebar/downloadmenu.R") #downloadmenu
source("gui/body/plots.R") #plots
source("gui/body/analysis.R") #plots

load('gui/obj/cell_ontology')

cell_ontology_names<-paste(cell_ont_full$id,cell_ont_full$name,sep = " ")

header <- dashboardHeader(
    titleWidth = 400,
    title = "Cellar"
)

sidebar <- dashboardSidebar(
    useShinyjs(),
    width = 400,
    sidebarMenu(
        datasetmenu,
        configmenu,
        selectionmenu,
        analysismenu,
        downloadmenu
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("gui/styles/style.css")),
    tags$script(src = "gui/sidebar/anchor.js"),
    plots,
    analysis
)

ui <- dashboardPage(header, sidebar, body, skin='green')
