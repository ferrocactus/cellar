pipe_sscluster <- function(pipe, new_labels, ssc_method, saved_clusters) {
    withProgress(message = "Running", value = 0, {
        incProgress(1 / 2, detail = paste("Constrained Clustering"))
        msg <- pipe()$run_step(step = "ssclu", ssclu_method = ssc_method,
                               ssclu_new_labels = new_labels,
                               saved_clusters = saved_clusters)
        return(msg)
    })
}

sscluster_run <- function(input, output, session, pipe, newLabels, replot) {
    observeEvent(input$ssclurun, {
        if (is.null(newLabels())) {
            showNotification("No updated plot found.")
            return()
        }

        msg <- pipe_sscluster(pipe, new_labels = newLabels(),
                              ssc_method = input$ssc_method,
                              saved_clusters = input$saved_clusters)

        if (msg != 'good') {
            showNotification(msg)
            return()
        }

        replot(1)
        #newLabels(pipe()$labels)
    })
}
