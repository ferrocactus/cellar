history <- function(id, label = "history") {
    ns = NS(id)
    fluidRow(
        div(class = "hist",
            conditionalPanel(
                'output.plot',
                ns = ns,
                actionButton(
                    ns("store_plot"), "Store Plot",
                    class = "histbtn"
                ),
                actionButton(
                    ns("delete_plot"), "Delete Plot",
                    class = "histbtn"
                )
            )
        )
    )
}
