info <- function(id, label="info") {
    ns <- NS(id)
    dropdownMenuOutput(ns("notifications"))
}
