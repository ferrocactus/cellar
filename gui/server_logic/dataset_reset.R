dataset_reset <- function(input, output, session, reset, setNames,
                          labelList, adata, fullreset, deGenes,
                          resubset, main_plot_val) {
    toListenReset <- reactive({
        list(reset(), fullreset())
    })

    observeEvent(toListenReset(), {
        if (reset() < 1 && fullreset() < 1) return()

        if (fullreset() > 0) {
            isolate(fullreset(0))
            isolate(reset(0))
            setNames(c())
            labelList(c())
            output$plot <- NULL
            main_plot_val(NULL)
        } else {
            isolate(reset(0))
        }

        updateTabsetPanel(session, "tabset", "Main Plot")

        output$DEtable = NULL
        output$GOtable = NULL
        output$KEGGtable = NULL
        output$CellTypetable = NULL
        output$UCellTypetable = NULL
        output$MSIGDBtable = NULL
        output$Diseasetable = NULL
        output$heatmap = NULL

        output$cell_names_outp = NULL
        output$clustering_info = NULL

        output$titleDE = NULL
        output$titleONTO = NULL
        output$titleKEGG = NULL
        output$titleMSIG = NULL
        output$titleCellType = NULL
        output$titleUCellType = NULL
        output$titleDisease = NULL
        output$titleheatmap = NULL

        updateSelectInput(
            session = session,
            inputId = "color",
            choices = c("Clusters"),
            selected = "Clusters")

        deGenes(c())

        req(adata())

        updateSliderInput(
            session = session,
            inputId = "mark_markers_n",
            min = 5, max = max(py_to_r(adata()$shape[1]), 10), value = 50
        )

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
                choices = c("Clusters","Uncertainty", names),
                selected = "Clusters")
        }

        # Update sets
        # anndata_has_key defined in r_helpers.py
        if (py_to_r(has_key(adata(), 'uns', 'cluster_names')) == TRUE) {
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

    })
}
