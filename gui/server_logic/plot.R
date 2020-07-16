all_symbols <- c(
    'star-triangle-down','circle','square',"diamond",
    "x",'triangle-up','triangle-down','hexagon','asterisk',
    'diamond-cross',"square-cross","circle-cross","circle-x",
    "star-square","star","star-triangle-up","star-square",
    "star-diamond","diamond-tall","diamond-wide","hourglass",
    'bowtie','pentagon','hexagram-dot','triangle-se','y-right',
    'hexagon2','octagon','triangle-nw','triangle-sw')

# Define the number of colors you want
plot <- function(input, output, session, replot, adata, selDataset,
                 setNames, setPts, plotHistory, curPlot, reset,
                 resubset, cellNamesTb, infoTb, reinfo, relabel) {
    # triggers when replot is set to 1
    observeEvent(replot(), {
        if (replot() < 1) return()
        isolate(replot(0))
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        relabel(relabel() + 1)
        reinfo(reinfo() + 1)

        labels <- as.factor(py_to_r(adata()$obs$labels$to_numpy()))
        label_names <- as.factor(py_to_r(get_label_names(adata())))
        x_emb_2d <- py_to_r(adata()$obsm[['x_emb_2d']])

        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            # Defaults
            title = "Clusters"
            showlegend = TRUE
            symbol = NULL
            symbols = NULL

            output$plot <- renderPlotly({
                text = ~paste("Label: ", label_names)
                if (input$color == 'Clusters') {
                    if (input$show_names == 'show_names')
                        color = paste0(labels, ": ", label_names)
                    else
                        color = labels
                } else {
                    gene_names = py_to_r(get_all_gene_names(adata()))
                    i = which(gene_names == input$color)[1]
                    if (is.null(i)) {
                        if (input$show_names == 'show_names')
                            color = paste0(labels, ": ", label_names)
                        else
                            color = labels
                    } else {
                        color = py_to_r(adata()$X)[, i]
                        title = input$color
                        showlegend = FALSE
                        symbol = ~labels
                        symbols = all_symbols
                    }
                }

                if (input$theme_mode == 'dark_mode') {
                    plot_bgcolor = 'rgb(44, 59, 65)'
                    paper_bgcolor = 'rgb(44, 59, 65)'
                    font = list(color = 'rgb(255, 255, 255)')
                    xaxis = list(gridcolor = 'rgb(85, 85, 85)',
                                 zerolinecolor = 'rgb(125, 125, 125)')
                    yaxis = list(gridcolor = 'rgb(85, 85, 85)',
                                 zerolinecolor = 'rgb(125, 125, 125)')
                } else {
                    plot_bgcolor = NULL
                    paper_bgcolor = NULL
                    font = NULL
                    xaxis = NULL
                    yaxis = NULL
                }

                p <- plotly_build(plot_ly(
                    x = x_emb_2d[, 1], y = x_emb_2d[, 2],
                    text = text,
                    color = color,
                    symbol = symbol,
                    symbols = symbols,  ## 30 shapes
                    key = as.character(1:length(labels)),
                    marker = list(size = input$dot_size),
                    type = 'scatter',
                    mode = 'markers',
                    height = input$plot_height
                ) %>% layout(dragmode = "lasso", showlegend = showlegend,
                             plot_bgcolor = plot_bgcolor,
                             paper_bgcolor = paper_bgcolor,
                             font = font,
                             xaxis = xaxis,
                             yaxis = yaxis,
                             title = title, margin = list(t = 50)))

                isolate(plotHistory(c(plotHistory(), list(p))))
                isolate(curPlot(length(plotHistory())))
                return(p)
            })
        })
    })

    # Listen to multiple events that will trigger replot
    toListenReplot <- reactive({
        list(
            input$plot_height,
            input$dot_size,
            input$color,
            input$show_names,
            input$theme_mode
        )
    })

    observeEvent(toListenReplot(), {
        replot(replot() + 1)
    })

    ###########################################################################
    # Plot Navigation
    observe({
        if (curPlot() < 1) return()

        output$plot <- renderPlotly({
            return(plotHistory()[[curPlot()]])
        })
        output$cell_names_outp <- renderUI({
            return(cellNamesTb()[[curPlot()]])
        })
        output$clustering_info <- renderUI({
            return(infoTb()[[curPlot()]])
        })
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
    ###########################################################################

    # Store selected cells
    observeEvent(input$store_lasso, {
        if (curPlot() == 0) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        if (curPlot() != length(plotHistory())) {
            showNotification("You can only select in the active plot.")
            return()
        }

        if (as.character(input$newsubset) %in% setNames()) {
            showNotification("Name already exists")
            return()
        }

        if (substr(as.character(input$newsubset), 1, 7) == "Cluster") {
            showNotification("Reserved name, please choose another.")
            return()
        }

        if (as.character(input$newsubset) == "") {
            showNotification("Please enter a name for the subset.")
            return()
        }

        d <- event_data("plotly_selected")
        if (is.null(d$key)) return()

        keys <- as.numeric(d$key)
        msg <- cellar$safe(store_subset,
            adata = adata(),
            name = input$newsubset,
            indices = keys,
            from_r = TRUE)

        if (msg != 'good') {
            showNotification(py_to_r(msg))
            return()
        }
        showNotification(paste(as.character(length(keys)), "cells stored"))
        resubset(resubset() + 1)
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
