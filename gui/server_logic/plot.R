# Define the number of colors you want
plot <- function(input, output, session, replot, pipe, selDataset,
                 setNames, setPts, plotHistory, curPlot, reset, relabel) {
    # triggers when replot is set to 1
    observeEvent(replot(), {
        if (replot() < 1) return()
        isolate(replot(0))
        if (pipe() == 0) return()
        if (!pipe()$has('labels')) return()


        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            label_names = pipe()$get_label_names()

            output$plot <- renderPlotly({
                if (input$color == 'Clusters') {
                    if (input$show_names == 'show_names')
                        color = paste0(
                            as.factor(pipe()$labels), ": ",
                            as.factor(label_names))
                    else
                        color = as.factor(pipe()$labels)
                    text = ~paste("Label: ", as.factor(label_names))
                    title = "Clusters"
                    p <- plotly_build(plot_ly(
                        x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                        text = text,
                        color = color,
                        key = as.character(1:length(pipe()$x_emb_2d[,1])),
                        marker = list(size = input$dot_size),
                        type = 'scatter',
                        mode = 'markers',
                    ) %>% layout(dragmode = "lasso", title = title,
                                 margin = list(t = 50)))



                } else {
                    i = which(pipe()$col_names == input$color)[1]
                    if (is.null(i)) {
                        if (input$show_names == 'show_names')
                            color = paste0(
                                as.factor(pipe()$labels), ": ",
                                as.factor(label_names))
                        else
                            color = as.factor(pipe()$labels)
                        title = "Clusters"
                    } else {
                        color = pipe()$x[, i]
                        title = input$color
                    }

                    text = ~paste("Label: ", as.factor(label_names))

                    p <- plotly_build(plot_ly(
                        x = pipe()$x_emb_2d[,1], y = pipe()$x_emb_2d[,2],
                        text = text,
                        color = color,
                        symbol = ~as.factor(pipe()$labels),
                        symbols = c('star-triangle-down','circle','square',"diamond",
                                    "x",'triangle-up','triangle-down','hexagon','asterisk','diamond-cross',
                                    "square-cross","circle-cross","circle-x","star-square","star","star-triangle-up",
                                    "star-square","star-diamond","diamond-tall","diamond-wide","hourglass",'bowtie',
                                    'pentagon','hexagram-dot','triangle-se','y-right','hexagon2','octagon','triangle-nw','triangle-sw'
                        ),  ## 30 shapes
                        key = as.character(1:length(pipe()$x_emb_2d[,1])),
                        marker = list(size = input$dot_size),
                        type = 'scatter',
                        mode = 'markers',
                    ) %>% layout(dragmode = "lasso", showlegend=FALSE,title = title,
                                 margin = list(t = 50)))

                }


                isolate(plotHistory(c(plotHistory(), list(p))))
                isolate(curPlot(length(plotHistory())))
                return(p)
            })
        })

    })

    observeEvent(input$prevplot, {
        if (curPlot() == 1) {
            showNotification("No more plots to show")
            return()
        }
        curPlot(curPlot() - 1)
        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
    })

    observeEvent(input$nextplot, {
        if (curPlot() == length(plotHistory())) {
            showNotification("This is the last plot")
            return()
        }
        curPlot(curPlot() + 1)
        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
    })

    # triggered when view gene expression is clicked
    observeEvent(input$color, {
        replot(replot() + 1)
    })

    # Store selected cells
    observeEvent(input$store_lasso, {
        if (curPlot() != length(plotHistory())) {
            showNotification("You can only select in the last plot.")
            return()
        }
        if (curPlot() == 0) return()
        if (!pipe()$has('labels')) return()

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
        if (pipe() == 0) return()
        if (!pipe()$has('labels')) return()

        idx <- which(setNames() == as.character(input$subset1_upd))
        if (idx == "None") return()

        keys <- setPts()[[idx]]
        # pass name and indices
        pipe()$update_labels(as.character(input$newlabels), keys - 1)
        reset(reset() + 1)
        replot(replot() + 1)
        relabel(relabel() + 1)
    })

    observeEvent(input$show_names, {
        replot(replot() + 1)
    })

    # Download plot
    output$download_plot <- downloadHandler(
        filename = function() {
            extension <- tolower(input$plot_download_format)
            paste0(tools::file_path_sans_ext(basename(selDataset())),
                   "_plot.", extension)
        },
        content = function(fname) {
            if (curPlot() == 0) return(NULL)
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
