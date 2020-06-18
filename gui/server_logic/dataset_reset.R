dataset_reset <- function(input, output, session, reset, setNames,
                          labelList, deButtons, deGenes, adata, fullreset,
                          curPlot, plotHistory, resubset) {
    observe({ if (reset() > 0 || fullreset() > 0) {
        if (fullreset() > 0) {
            setNames(c("None"))
            labelList(c())
            curPlot(0)
            plotHistory(c())
            output$plot <- NULL
        } else {
            resubset(resubset() + 1)
        }

        isolate(reset(0))
        isolate(fullreset(0))

        updateSelectInput(
            session = session,
            inputId = "color",
            choices = c("Clusters", as.character(as.character(adata()$var_names))),
            selected = "Clusters")

        # Update sets
        # anndata_has_key defined in r_helpers.py
        if (has_key(adata(), 'obs', 'labels') == TRUE) {
            resubset(resubset() + 1)

            if (input$tissue == 'clusters') {
                updateSelectInput(
                    session,
                    "newlabels",
                    choices = n_clusters
                )
            }
        } else {
            if (input$tissue == 'clusters') {
                updateSelectInput(
                    session,
                    "newlabels",
                    choices = c()
                )
            }
        }

        # Clear all analysis tabs
        if (length(deButtons()) > 0)
            for (i in 1:length(deButtons()))
                deButtons()[[i]]$destroy()
        deButtons(c())
        deGenes(c())

        output$DEbuttons = NULL
        output$DEtable = NULL
        output$GOtable = NULL
        output$KEGGtable = NULL
        output$Markerstable = NULL
        output$Markerstableuser = NULL
        output$MSIGDBtable = NULL

    }})
}
