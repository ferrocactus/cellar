source_python("gui/server_logic/read_onto.py")

update_label <- function(input, output, session, pipe, labelList) {
    dic = get_dic()

    observeEvent(input$tissue, {
        tissue = as.character(input$tissue)
        if (tissue == "") return()

        if (tissue == 'clusters') {
            if (pipe() == 0) return()

            n_clusters = pipe()$n_clusters
            updateSelectInput(
                session,
                "newlabels",
                choices = n_clusters
            )
        } else if (tissue == 'user defined') {
            if (length(labelList()) == 0) return()

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
        labelList(union(labelList(), input$newlabelbox))
        if (input$tissue == 'user defined') {
            updateSelectInput(
                session,
                "newlabels",
                choices = labelList()
            )
        }
        showNotification("Label added")
    })
}

