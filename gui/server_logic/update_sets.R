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
            "subset1_upd",
            "Choose Subset",
            choices = setNames()
    )})

    # observe({
    #     updateTextInput(
    #         session,
    #         "cell_subset_download",
    #         "Input Subset IDs",
    #         value='0'
    #         #value=strsplit(setNames()[2],'_')[[1]][[2]]
    # )})
}
