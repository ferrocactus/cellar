dataset_reset <- function(input, output, session, reset, setNames,
                          labelList, deGenes, adata, fullreset,
                          curPlot, plotHistory, resubset) {
    toListenReset <- reactive({
        list(reset(), fullreset())
    })

    observeEvent(toListenReset(), {
        if (reset() < 1 && fullreset() < 1) return()

        if (fullreset() > 0) {
            isolate(fullreset(0))
            isolate(reset(0))
            setNames(c("None"))
            labelList(c())
            curPlot(0)
            plotHistory(c())
            output$plot <- NULL
        } else {
            isolate(reset(0))
            resubset(resubset() + 1)
        }

        updateSelectInput(
            session = session,
            inputId = "color",
            choices = c("Clusters", as.character(as.character(adata()$var_names))),
            selected = "Clusters")

        # Update sets
        # anndata_has_key defined in r_helpers.py
        if (has_key(adata(), 'uns', 'cluster_names') == TRUE) {
            resubset(resubset() + 1)

            clusters = py_to_r(get_cluster_label_list(adata()))
            if (input$tissue == 'clusters') {
                updateSelectInput(
                    session,
                    "newlabels",
                    choices = clusters
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

        deGenes(c())

        output$DEtable = NULL
        output$GOtable = NULL
        output$KEGGtable = NULL
        output$Markerstable = NULL
        output$Markerstableuser = NULL
        output$MSIGDBtable = NULL

    })
}
