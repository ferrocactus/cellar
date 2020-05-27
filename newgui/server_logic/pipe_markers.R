# Needed by markers_run
pipe_markers <- function(pipe, set1, set2, mark_method, mark_alpha,
                         mark_markers_n, mark_correction, con_method,
                         con_convention, con_path) {
    withProgress(message = "DE", value = 0, {
       n <- 1
       incProgress(1/n, detail = paste("Step: Finding DE genes"))
       pipe()$get_markers_subset(indices1 = set1, indices2 = set2,
                                 method = mark_method, alpha = mark_alpha,
                                 markers_n = mark_markers_n,
                                 correction = mark_correction,
                                 con_method = con_method,
                                 con_convention = con_convention,
                                 con_path = con_path)
    })
}

markers_run <- function(input, output, session, pipe, remark, setNames, setPts) {
    observeEvent(input$getdegenes, {
        s1 = as.character(input$subset1)
        s2 = as.character(input$subset2)

        if (s1 == 'None' && s2 == 'None') {
            showNotification("Please select a subset to analyze.")
            return()
        }

        if (s1 == 'None') {
            set1 = setPts()[[which(setNames() == s2)]]
            set2 = NULL
        } else if (s2 == 'None') {
            set1 = setPts()[[which(setNames() == s1)]]
            set2 = NULL
        } else {
            set1 = setPts()[[which(setNames() == s1)]]
            set2 = setPts()[[which(setNames() == s2)]]
        }

        # Subtract by 1 due to R lists starting at 1
        pipe_markers(pipe, set1 = (set1 - 1), set2 = (set2 - 1),
                    mark_method = 'TTest',
                    mark_alpha = input$mark_alpha,
                    mark_markers_n = input$mark_markers_n,
                    mark_correction = input$mark_correction,
                    con_method = 'Converter',
                    con_convention = 'id-to-name',
                    con_path = 'markers/gene_id_name.csv')

        remark(remark() + 1) # Notify that markers have changed
        #TODO Update tabsetpanel to switch to DE
        #ns <- session$ns
        #updateTabsetPanel(session, ns("switcher"), selected = ns("DE"))
    })
}
