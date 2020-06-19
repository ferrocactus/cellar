dataset_reset <- function(input, output, session, reset, setNames,
                          labelList, adata, fullreset, deGenes,
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
        }

        deGenes(c())

        if (is_active(adata()) == FALSE) return()

        if (has_key(adata(), 'var', 'parsed_names') == FALSE) {
            updateSelectInput(
                session = session,
                inputId = "color",
                choices = c())
        } else {
            names = py_to_r(get_all_gene_names(adata()))
            updateSelectInput(
                session = session,
                inputId = "color",
                choices = c("Clusters", names),
                selected = "Clusters")
        }

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

        output$DEtable = NULL
        output$GOtable = NULL
        output$KEGGtable = NULL
        output$CellTypetable = NULL
        output$UCellTypetable = NULL
        output$MSIGDBtable = NULL
        output$Diseasetable = NULL
        output$heatmap = NULL
    })
}