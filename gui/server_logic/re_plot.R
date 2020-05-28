re_plot <- function(input, output, session, replot, pipe, plotObj, selDataset) {
    observe({ if (replot() > 0) {
        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            output$plot <- renderPlotly({
                # TODO fix colors and key
                if (input$color == "cluster") {
                    plotcols = as.factor(pipe()$labels)
                } else {
                    plotcols = pipe()$labels
                }

                plotObj(plot_ly(
                    x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                    text = ~paste("Label: ", as.factor(pipe()$labels)),
                    color = as.factor(pipe()$labels),
                    key = plotcols,
                    type = 'scatter',
                    mode = 'markers'
                ) %>% layout(dragmode = "lasso",
                            title = paste("Value of ", input$color, sep="")))
                print(plotObj())
                return(plotObj())
            })
        })
        replot(0)
    }})


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
