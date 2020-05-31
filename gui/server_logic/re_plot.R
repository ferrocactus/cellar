re_plot <- function(input, output, session, replot, pipe, plotObj, selDataset) {
    observe({ if (replot() > 0) {
        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            output$plot <- renderPlotly({
                if (input$color == 'Clusters') {
                    color = as.factor(pipe()$labels)
                    title = "Clusters"
                } else {
                    i = which(pipe()$col_ids == input$color)[1]
                    if (is.null(i)) {
                        color = as.factor(pipe()$labels)
                    } else {
                        color = pipe()$x[, i]
                        title = input$color
                    }
                }

                plotObj(plot_ly(
                    x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                    text = ~paste("Label: ", as.factor(pipe()$labels)),
                    color = color,
                    key = as.character(1:length(pipe()$x_emb_2d[,1])),
                    type = 'scatter',
                    mode = 'markers'
                ) %>% layout(dragmode = "lasso",
                            title = title))
                return(plotObj())
            })
        })
        replot(0)
    }})

    observeEvent(input$color, suspended=TRUE, {
        replot(1)
    })

    # Download plot
    output$download_plot <- downloadHandler(
        filename = function() {
            extension <- tolower(input$plot_download_format)
            paste0(tools::file_path_sans_ext(basename(selDataset())),
                   "_plot.", extension)
        },
        content = function(fname) {
            extension <- tolower(input$plot_download_format)
            if (extension == 'html') {
                htmlwidgets::saveWidget(as_widget(plotObj()), fname,
                                    selfcontained = TRUE)
            } else {
                withProgress(message = "Saving plot", value = 0, {
                    incProgress(1/2, detail = paste("Step: Rendering"))
                    withr::with_dir(dirname(fname), orca(plotObj(),
                            basename(fname), format = extension))
                })
            }
        }
      )

}
