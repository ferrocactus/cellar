add_label <- function(input, output, session, labelList) {
    observe({
        updateSelectInput(
            session,
            "newlabels",
            label = paste("Select input label", length(x)),
            choices = labelList())
    })

    observeEvent(input$labeladd, {
        labelList(union(labelList(), input$newlabelbox))
    })
}
