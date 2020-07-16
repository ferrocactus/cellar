cell_names <- function(id, label='cell_names') {
    ns = NS(id)
    div(
        class = "cell_names_div",
        list(
            actionButton(
                ns("collapse_cell_names"),
                "View additional info",
                class = "collapsebtn"
            ),
            splitLayout(
                htmlOutput(ns("clustering_info")),
                htmlOutput(ns("cell_names_outp"))
            )
        )
    )
}
