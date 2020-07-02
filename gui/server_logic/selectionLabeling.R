

selectionLabeling <- function(input, output, session, adata,
                              labelList, setNames, resubset,
                              reset, replot, relabel) {

    dic = py_to_r(get_dic())

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
        if (input$newlabelbox == '') return()
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
        if (input$newlabels == 'None') return()
        if (input$subset1_upd == 'None') return()

        msg <- cellar$safe(update_subset_label,
            adata = adata(),
            subset_name = input$subset1_upd,
            name = input$newlabels)

        if (msg != 'good') {
            showNotification(py_to_r(msg))
            return()
        }
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
    
    # observe({
    #     print("rendering")
    #     output$subset1_upd<-renderUI({
    #     selectInput(
    #         ns("subset1_upd"),
    #         "Choose Subset",
    #         choices = c("None", setNames())
    #     ) 
        # %>%
        #     shinyInput_label_embed(
        #         shiny::icon("info-circle") %>%
        #             bs_embed_tooltip(
        #                 "Select the subset to be updated",
        #                 placement= "bottom"
        #             )
    #     #     )
    #     })
    # })
    
    observe({
        updateSelectInput(
            session,
            "subset1_upd",
            "Choose Subset",
            choices = c("None", setNames())
        )
    })
}