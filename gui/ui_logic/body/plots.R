plots <- function(id, label='plots') {
ns = NS(id)
tabsetPanel(
    type = "tabs",
    id = ns("tabset"),
    tabPanel(
        "Main Plot",
        uiOutput(ns("plots")),

        conditionalPanel(
            "output.plot",
            ns = ns,
            div(
                class = "cell_names_div",
                list(
                    actionButton(
                        ns("split_plot"),
                        "Split"
                    ),
                    actionButton(
                        ns("collapse_cell_names"),
                        "View additional info",
                        class = "collapsebtn"
                    ),
                    splitLayout(
                        htmlOutput(ns("clustering_info")),
                        htmlOutput(ns("cell_names_outp"))
            )))
        )
    )
)}
