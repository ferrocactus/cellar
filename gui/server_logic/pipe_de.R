# Needed by markers_run
pipe_de <- function(pipe, set1, set2, mark_method, mark_alpha,
                    mark_markers_n, mark_correction, con_method,
                    con_convention, con_path) {
    withProgress(message = "DE", value = 0, {
        n <- 2
        incProgress(1 / n, detail = paste("Step: Finding DE genes"))
        msg <- pipe()$run_step(step = 'de_subset', subset1 = set1, subset2 = set2,
                               de_method = mark_method, de_alpha = mark_alpha,
                               de_n_genes = mark_markers_n,
                               de_correction = mark_correction)
        if (msg != 'good') return(msg)
        incProgress(1 / n, detail = paste("Step: Converting names"))
        msg <- pipe()$run_step(step = 'con', con_method = con_method,
                               con_convention = con_convention,
                               con_path = con_path)
        return(msg)
    })
}

de_run <- function(input, output, session, pipe, remark, setNames, setPts) {
    observeEvent(input$getdegenes, {
        if (pipe() == 0) return()

        s1 = as.character(input$subset1)
        s2 = as.character(input$subset2)

        if (s1 == 'None' && s2 == 'None') {
            showNotification("Please select a subset to analyze.")
            return()
        }

        if (s1 == s2) {
            showNotification("The selected subsets should be different.")
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
        msg <- pipe_de(pipe, set1 = (set1 - 1), set2 = (set2 - 1),
                    mark_method = 'TTest',
                    mark_alpha = input$mark_alpha,
                    mark_markers_n = input$mark_markers_n,
                    mark_correction = input$mark_correction,
                    con_method = 'Converter',
                    con_convention = 'id-to-name',
                    con_path = 'markers/gene_id_name.csv')

        if (msg != 'good') {
            showNotification(msg)
        } else {
            remark(remark() + 1) # Notify that markers have changed
        }
    })
}
