dataset_reset <- function(input, output, session, reset, setNames, setPts,
                          labelList, deButtons, deGenes, adata, fullreset,
                          curPlot, plotHistory) {
    observe({ if (reset() > 0 || fullreset() > 0) {
        if (fullreset() > 0) {
            setNames(c("None"))
            setPts(c(NA))
            labelList(c())
            curPlot(0)
            plotHistory(c())
            output$plot <- NULL
        } else {
            setPts(setPts()[which(substr(setNames(), 1, 7) != 'Cluster')])
            setNames(setNames()[which(substr(setNames(), 1, 7) != 'Cluster')])
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
            n_clusters = adata()$uns$cluster_info$n_clusters
            for (i in n_clusters) {
                setNames(c(setNames(),
                           (paste("Cluster_", as.character(i), sep = ""))))
                setPts(c(setPts(), list(which(adata()$obs$labels == i))))
            }

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
