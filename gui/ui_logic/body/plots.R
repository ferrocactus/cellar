library(shiny)

plots <- function(id, label='plots') {
ns = NS(id)
tabsetPanel(
    type = "tabs",
    id = "tabset",
    tabPanel(
        "Main Plot",
        h3(textOutput("caption")),
        plotlyOutput(ns("plot"), height="600px"),
    )
)}
