de <- function(input, output, session, pipe) {
    if (as.character(input$subset2) == "None"
        && as.character(input$subset1) == "None") {
        showNotification("Please select at least a subset to run the analysis.")
        return()
    }

    selecteddat = NULL
    idx1 = which(set_name == as.character(input$subset1))
    idx2 = which(set_name == as.character(input$subset2))

    assign("s1", sets[[idx1]], envir = env)
    assign("s2", sets[[idx2]], envir = env)

    cell_count=length(s1)
    cell_count2=length(s2)
}
