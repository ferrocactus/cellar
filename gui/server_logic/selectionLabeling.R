source_python("gui/server_logic/read_onto.py")

selectionLabeling <- function(input, output, session, adata,
                              labelList, setNames, resubset,
                              reset, replot, relabel) {
    dic = get_dic()

    observeEvent(input$tissue, {
        tissue = as.character(input$tissue)
        if (tissue == "") return()

        if (tissue == 'clusters') {
            if (is_active(adata()) == FALSE) return()
            if (!py_has_attr(adata()$obs, 'labels')) return()

            clusters = py_to_r(get_unique_labels(adata()))
            updateSelectInput(
                session,
                "newlabels",
                choices = clusters
            )
        } else if (tissue == 'user defined') {
            updateSelectInput(
                session,
                "newlabels",
                choices = labelList()
            )
        } else if (tissue == 'all') {
            choices = c()
            for (tissue in names(dic)) {
                choices = c(choices, dic[tissue])
            }
            updateSelectInput(
                session,
                "newlabels",
                choices = choices
            )
        } else {
            updateSelectInput(
                session,
                "newlabels",
                choices = dic[tissue]
            )
        }
    })

    # Observe adding new labels
    observeEvent(input$labeladd, {
        if (is.null(input$newlabelbox)) return()
        labelList(union(labelList(), input$newlabelbox))
        showNotification("Label added")
    })

    observeEvent(labelList(), {
        if (input$tissue == 'user defined') {
            updateSelectInput(
                session,
                "newlabels",
                choices = labelList()
            )
        }
    })

    # Observe change of subsets
    observeEvent(resubset(), {
        if (resubset() < 1) return()
        isolate(resubset(0))

        setNames(c())

        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        subsets <- py_to_r(get_subsets(adata()))
        for (subset in subsets)
            setNames(c(setNames(), as.character(subset)))
    })

    # Update the label of a subset
    observeEvent(input$labelupd, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        update_subset_label(adata(), input$subset1_upd, input$newlabels)
        reset(reset() + 1)
        replot(replot() + 1)
        relabel(relabel() + 1)
    })

    # Observe changes of sets
    observe({
        updateSelectInput(
            session,
            "subset1",
            "Choose Subset 1",
            choices = c("None", setNames())
    )})

    observe({
        updateSelectInput(
            session,
            "subset2",
            "Choose Subset 2",
            choices = c("None", setNames())
    )})

    observe({
        updateSelectInput(
            session,
            "subset1_upd",
            "Choose Subset",
            choices = setNames()
    )})
}