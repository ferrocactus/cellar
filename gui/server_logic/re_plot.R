re_plot <- function(input, output, session, replot, pipe, plotObj, selDataset,
                    setNames, setPts) {
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

                ns <- session$ns
                plotObj(plot_ly(
                    x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                    text = ~paste("Label: ", as.factor(pipe()$labels)),
                    color = color,
                    key = as.character(1:length(pipe()$x_emb_2d[,1])),
                    type = 'scatter',
                    mode = 'markers',
                    source = ns('M')
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

    observeEvent(input$store_lasso, {
        if (as.character(input$newsubset) %in% setNames()) {
          showNotification("Name already exists")
          return()
        }

        if (substr(as.character(input$newsubset),1,7) == "cluster") {
          showNotification("Reserved name, please choose another.")
          return()
        }

        ns <- session$ns
        d <- event_data("plotly_selected", source='M')
        keys <- as.numeric(d$key)
        cell_count <- length(d$key)
        print(input$newsubset)
        print(d$key)

        #showNotification(as.character(input$newsubset),duration=NULL)
        if (identical(d$key, NULL) == TRUE) {
          return()
        }

        if (as.character(input$newsubset) != "") {
          setPts(c(setPts(), list(keys)))
          setNames(c(setNames(), as.character(input$newsubset)))
        } else {
          setPts(c(setPts(), list(keys)))
          setNames(c(setNames(),
                    paste("Newset", as.character(length(setNames())), sep="")))
        }

        showNotification(paste(as.character(cell_count), "cells stored"))
    })
}
