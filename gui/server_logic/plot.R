symbols <- c(
    'star-triangle-down','circle','square',"diamond",
    "x",'triangle-up','triangle-down','hexagon','asterisk',
    'diamond-cross',"square-cross","circle-cross","circle-x",
    "star-square","star","star-triangle-up","star-square",
    "star-diamond","diamond-tall","diamond-wide","hourglass",
    'bowtie','pentagon','hexagram-dot','triangle-se','y-right',
    'hexagon2','octagon','triangle-nw','triangle-sw')

# Define the number of colors you want
plot <- function(input, output, session, replot, adata, selDataset,
                 setNames, setPts, plotHistory, curPlot, reset, relabel) {
    # triggers when replot is set to 1
    observeEvent(replot(), {
        if (replot() < 1) return()
        isolate(replot(0))
        if (is_active(adata()) == FALSE) return()
        if (has_key(adata(), 'obs', 'labels') == FALSE) return()

        labels <- as.factor(py_to_r(adata()$obs$labels$to_numpy()))
        x_emb_2d <- py_to_r(adata()$obsm[['x_emb_2d']])

        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            label_names = as.factor(py_to_r(get_label_names(adata())))

            output$plot <- renderPlotly({
                if (input$color == 'Clusters') {
                    if (input$show_names == 'show_names')
                        color = paste0(
                            labels, ": ",
                            label_names)
                    else
                        color = labels
                    text = ~paste("Label: ", label_names)
                    title = "Clusters"

                    p <- plotly_build(plot_ly(
                        x = x_emb_2d[, 1], y = x_emb_2d[, 2],
                        text = text,
                        color = color,
                        key = as.character(1:length(labels)),
                        marker = list(size = input$dot_size),
                        type = 'scatter',
                        mode = 'markers',
                        height = input$plot_height
                    ) %>% layout(dragmode = "lasso", title = title,
                                 margin = list(t = 50)))
                } else {
                    i = which(py_to_r(adata()$var_names$to_numpy()) == input$color)[1]
                    if (is.null(i)) {
                        if (input$show_names == 'show_names')
                            color = paste0(
                                labels, ": ",
                                label_names)
                        else
                            color = labels
                        title = "Clusters"
                    } else {
                        color = py_to_r(adata()$X)[, i]
                        title = input$color
                    }

                    text = ~paste("Label: ", label_names)

                    p <- plotly_build(plot_ly(
                        x = x_emb_2d[, 1], y = x_emb_2d[, 2],
                        text = text,
                        color = color,
                        symbol = ~labels,
                        symbols = symbols,  ## 30 shapes
                        key = as.character(1:length(labels)),
                        marker = list(size = input$dot_size),
                        type = 'scatter',
                        mode = 'markers',
                        height = input$plot_height
                    ) %>% layout(dragmode = "lasso", showlegend=FALSE,title = title,
                                 margin = list(t = 50)))
                }

                isolate(plotHistory(c(plotHistory(), list(p))))
                isolate(curPlot(length(plotHistory())))
                return(p)
            })
        })

    })

    observe({
        if (curPlot() < 1) return()

        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
    })

    observeEvent(input$plot_height, {
        replot(replot() + 1)
    })

    observeEvent(input$dot_size, {
        replot(replot() + 1)
    })

    observeEvent(input$prevplot, {
        if (curPlot() < 1) return()

        if (curPlot() == 1) {
            showNotification("No more plots to show")
            return()
        }
        curPlot(curPlot() - 1)
    })

    observeEvent(input$nextplot, {
        if (curPlot() < 1) return()

        if (curPlot() == length(plotHistory())) {
            showNotification("This is the last plot")
            return()
        }
        curPlot(curPlot() + 1)
    })

    observeEvent(input$firstplot, {
        if (curPlot() < 1) return()
        curPlot(1)
    })

    observeEvent(input$lastplot, {
        if (curPlot() < 1) return()
        curPlot(length(plotHistory()))
    })

    # triggered when view gene expression is clicked
    observeEvent(input$color, {
        replot(replot() + 1)
    })

    # Store selected cells
    observeEvent(input$store_lasso, {
        if (curPlot() != length(plotHistory())) {
            showNotification("You can only select in the active plot.")
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
        if (idx == 1) return()

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
