history <- function(id, label = "history") {
    ns = NS(id)
    fluidRow(
        div(class = "hist",
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
}
