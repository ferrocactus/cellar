cell_names <- function(id, label='cell_names') {
    ns = NS(id)
    div(
        class = "cell_names_div",
        list(
            actionButton(
                ns("collapse_cell_names"),
                "\u2199 View cluster names",
                class = "collapsebtn"
            ),
            tableOutput(ns("cell_names_outp"))
        )
    )
}
