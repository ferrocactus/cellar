re_label <- function(input, output, session, relabel, pipe) {
    observeEvent(input$collapse_cell_names, {
        ns <- session$ns
        shinyjs::toggle("cell_names_outp")
    })
    observe({
        if (relabel() < 1) return()
        output$cell_names_outp <- renderTable(width = "100%", {
            d <- pipe()$get_cluster_names()
            df <- data.frame(d$'labels', d$'names')
            colnames(df) <- c("Cluster ID", "Name")
            isolate(relabel(0))
            return(df)
        })
    })
}
