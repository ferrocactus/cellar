library(shiny)
library(shinydashboard)
library(plotly)
library(shinyjs)

source("gui/sidebar/datasetmenu.R") #datasetmenu
source("gui/sidebar/mainmenu.R") #mainmenu
source("gui/sidebar/configmenu.R") #configmenu
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
        mainmenu,
        configmenu
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("gui/styles/style.css")),
    # Don't scroll to top after clicking anything
    #tags$script(HTML(
    #    "
    #    $(document).on('click', function(e) {
    #        return false;
    #        e.preventDefault();
    #    });
    #    "
    #)),
    plots,
    analysis
)

ui <- dashboardPage(header, sidebar, body, skin='green')
