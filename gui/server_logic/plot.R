# Define the number of colors you want
plot <- function(input, output, session, replot, pipe, selDataset,
                 setNames, setPts, plotHistory, curPlot, reset, relabel) {
    # triggers when replot is set to 1
    observe({
        if (is.null(replot())) return()
        if (!pipe()$has('labels')) return()

        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            output$plot <- renderPlotly({
                if (input$color == 'Clusters') {
                    if (input$show_names == 'show_names')
                        color = paste0(
                                    as.factor(pipe()$labels), ": ",
                                    as.factor(pipe()$label_names))
                    else
                        color = as.factor(pipe()$labels)
                    title = "Clusters"
                } else {
                    i = which(pipe()$col_names == input$color)[1]
                    if (is.null(i)) {
                        if (input$show_names == 'show_names')
                            color = paste0(
                                        as.factor(pipe()$labels), ": ",
                                        as.factor(pipe()$label_names))
                        else
                            color = as.factor(pipe()$labels)
                        title = "Clusters"
                    } else {
                        color = pipe()$x[, i]
                        title = input$color
                    }
                }

                text = ~paste("Label: ", as.factor(pipe()$label_names))
                vals <- schema(F)$traces$scatter$attributes$marker$symbol$values
                vals <- grep("-", vals, value = T)
                p <- plotly_build(plot_ly(
                    x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                    text = text,
                    color = color,
                    symbol = ~as.factor(pipe()$labels), symbols = vals,#c('circle','x','o'),
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

        if (substr(as.character(input$newsubset),1,7) == "Cluster") {
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
        idx <- which(setNames() == as.character(input$subset1_upd))
        keys <- setPts()[[idx]]
        # pass name and indices
        pipe()$update_labels(as.character(input$newlabels), keys - 1)
        reset(1)
        replot(1)
        relabel(1)
    })

    observeEvent(input$show_names, suspended=TRUE, {
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
