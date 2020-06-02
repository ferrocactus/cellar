library(RColorBrewer)
# Define the number of colors you want
plot <- function(input, output, session, replot, pipe, selDataset,
                 setNames, setPts, newLabels, plotHistory, curPlot) {
    # triggers when replot is set to 1
    observe({
        req(replot())
        if (!pipe()$has('labels')) return()

        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            output$plot <- renderPlotly({
                if (input$color == 'Clusters') {
                    color = as.factor(pipe()$labels)
                    title = "Clusters"
                } else {
                    i = which(pipe()$col_names == input$color)[1]
                    if (is.null(i)) {
                        color = as.factor(pipe()$labels)
                        title = "Clusters"
                    } else {
                        color = pipe()$x[, i]
                        title = input$color
                    }
                }

                p <- plotly_build(plot_ly(
                    x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                    text = ~paste("Label: ", as.factor(pipe()$labels)),
                    color = color,
                    key = as.character(1:length(pipe()$x_emb_2d[,1])),
                    marker = list(size = input$dot_size),
                    type = 'scatter',
                    mode = 'markers',
                ) %>% layout(dragmode = "lasso", title = title,
                            margin = list(t = 50)))
                isolate(plotHistory(c(plotHistory(), list(p))))
                isolate(curPlot(length(plotHistory())))
                return(p)
            })
        })
        isolate(replot(NULL))
    })

    observeEvent(input$prevplot, {
        if (curPlot() == 1) return()
        curPlot(curPlot() - 1)
        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
    })

    observeEvent(input$nextplot, {
        if (curPlot() == length(plotHistory())) return()
        curPlot(curPlot() + 1)
        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
    })

    # triggered when view gene expression is clicked
    observeEvent(input$color, suspended=TRUE, {
        replot(1)
    })

    # Store selected cells
    observeEvent(input$store_lasso, {
        if (as.character(input$newsubset) %in% setNames()) {
          showNotification("Name already exists")
          return()
        }

        if (substr(as.character(input$newsubset),1,7) == "cluster") {
          showNotification("Reserved name, please choose another.")
          return()
        }

        d <- event_data("plotly_selected")
        keys <- as.numeric(d$key)
        cell_count <- length(keys)
        keys <- list(keys)

        if (is.null(d$key))
            return()

        if (as.character(input$newsubset) != "") {
            setPts(c(setPts(), keys))
            setNames(c(setNames(), as.character(input$newsubset)))
        } else {
            setPts(c(setPts(), keys))
            setNames(c(setNames(),
                paste("Newset", as.character(length(setNames())), sep="")))
        }

        showNotification(paste(as.character(cell_count), "cells stored"))
    })

    observeEvent(input$labelupd, {
        if (is.null(newLabels())) isolate(newLabels(pipe()$labels))
        idx <- which(setNames() == as.character(input$subset1_upd))
        keys <- setPts()[[idx]]
        lbs <- newLabels()
        lbs[keys] <- as.character(input$newlabels)
        isolate(newLabels(lbs))
        showNotification("Updating labels")

        if (input$color == 'Clusters')
            title = "Clusters"
        else
            title = input$labelupd

        output$plot <- renderPlotly({
            p <- plotly_build(plot_ly(
                x = pipe()$x_emb_2d[, 1], y = pipe()$x_emb_2d[, 2],
                text = ~paste("label: ", as.factor(newLabels())),
                color = as.factor(newLabels()),
                key = as.character(1:length(pipe()$x_emb_2d[,1])),
                type = 'scatter',
                mode = 'markers',
            ) %>% layout(dragmode = "lasso", title = title,
                        margin = list(t = 50)))
            isolate(plotHistory(c(plotHistory(), list(p))))
            isolate(curPlot(length(plotHistory())))
            return(p)
        })
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
                htmlwidgets::saveWidget(
                        as_widget(plotHistory()[[curPlot()]]), fname,
                        selfcontained = TRUE)
            } else {
                withProgress(message = "Saving plot", value = 0, {
                    incProgress(1/2, detail = paste("Step: Rendering"))
                    withr::with_dir(dirname(fname),
                            orca(plotHistory()[[curPlot()]],
                            basename(fname), format = extension))
                })
            }
        }
    )
}
