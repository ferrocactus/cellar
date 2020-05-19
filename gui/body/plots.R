library(shiny)

plots <- tabsetPanel(
    type = "tabs",
    id = "tabset",
    tabPanel(
        "Main Plot",
        h3(textOutput("caption")),
        plotlyOutput("plot", height="550px"),
    ),
    tabPanel(
        "Updated Plot",
        verbatimTextOutput("brush"),
        plotlyOutput("Plot2", height="550px")
    )
)
