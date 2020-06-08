history <- function(id, label = "history") {
    ns = NS(id)
    fluidRow(
        div(class = "hist",
            actionButton(
                ns("firstplot"), icon("angle-double-left"),
                class = "histbtn prevbtn"
            ),
            actionButton(
                ns("prevplot"), icon("angle-left"),
                class = "histbtn prevbtn"
            ),
            actionButton(
                ns("nextplot"), icon("angle-right"),
                class = "histbtn nextbtn"
            ),
            actionButton(
                ns("lastplot"), icon("angle-double-right"),
                class = "histbtn nextbtn"
            )
        )
    )
}
