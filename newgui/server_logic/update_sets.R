update_sets <- function(input, output, session, setNames) {
    observe({
        updateSelectInput(
            session,
            "subset1",
            "Choose Subset 1",
            choices = setNames()
    )})

    observe({
        updateSelectInput(
            session,
            "subset2",
            "Choose Subset 2",
            choices = setNames()
    )})

    observe({
        updateSelectInput(
            session,
            "cell_subset_download",
            "Download labels",
            choices = c("All", setNames())
    )})
}
