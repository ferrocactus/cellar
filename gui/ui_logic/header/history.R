history <- function(id, label = "history") {
    ns = NS(id)
    fluidRow(
        div(class = "hist",
            actionButton(
                ns("prevplot"), "Prev",
                class = "histbtn prevbtn"
            ),
            actionButton(
                ns("nextplot"), "Next",
                class = "histbtn nextbtn"
            )
        )
    )
}
