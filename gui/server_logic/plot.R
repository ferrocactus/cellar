source("gui/server_logic/plot_options.R")

# Define the number of colors you want
plot <- function(input, output, session, replot, adata, selDataset,
                 setNames, setPts, plotHistory, curPlot, reset,
                 resubset, cellNamesTb, infoTb, reinfo, relabel) {
    no_append = reactiveVal(0)
    navigate = reactiveVal(0)

    # triggers when replot is set to 1
    observeEvent(replot(), {
        if (replot() < 1) return()
        isolate(replot(0))
        no_app <- isolate(no_append())
        isolate(no_append(0))

        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        if (no_app == 0) {
            relabel(relabel() + 1)
            reinfo(reinfo() + 1)
        }

        withProgress(message = "Making plot", value = 0, {
            incProgress(1, detail = paste("Step: Rendering plot"))
            labels <- as.factor(py_to_r(adata()$obs$labels$to_numpy()))
            label_names <- as.factor(py_to_r(get_label_names(adata())))
            x_emb_2d <- py_to_r(adata()$obsm[['x_emb_2d']])

            title = "Clusters"
            showlegend = TRUE
            symbol = NULL
            symbols = NULL

            text = ~paste("Label: ", label_names)
            if (input$show_names == 'show_names' && input$color == 'Clusters')
                color = paste0(labels, ": ", label_names)
            else if (input$color == 'Clusters')
                color = labels
            else if (isolate(input$color) == 'Uncertainty'){
                if (anyNA(as.integer(isolate(input$n_neighbors)))==TRUE){
                    n_neighbors=as.integer(sqrt(length(labels)))
                } else if (as.integer(isolate(input$n_neighbors)) > as.integer(length(labels)/2)){
                    n_neighbors=as.integer(sqrt(length(labels)))
                } else{
                    n_neighbors=as.integer(isolate(input$n_neighbors))
                }

                showNotification("Calculating Uncertainty")
                withProgress(message = "Please Wait", value = 0, {
                    incProgress(1 / 3, detail = "Data processing...")
                    cellar$safe(cellar$get_neighbors,
                                x=adata(),
                                n_neighbors=n_neighbors)
                    incProgress(1 / 3, detail = "Calculating uncertainty...")
                    cellar$safe(cellar$uncertainty,
                                x=adata(),
                                n_neighbors=n_neighbors)
                })
                showNotification("Finished")

                color = as.numeric(py_to_r(adata()$obs['uncertainty']))
                title = "Uncertainty"
                showlegend = FALSE
                symbol = ~labels
                symbols = all_symbols
                cell_ids = as.character(py_to_r(adata()$obs$index$to_numpy()))
                certainty = py_to_r(adata()$uns['uncertainty_text'])
                text = ~paste0(replicate(length(label_names),"Label: "),
                            label_names,replicate(length(label_names),", Uncertainty: "),
                            as.character(py_to_r(adata()$obs['uncertainty'])),'\n',as.character(certainty))
            } else {
                gene_names = py_to_r(get_all_gene_names(adata()))
                i = which(gene_names == isolate(input$color))[1]
                if (!is.null(i)) {
                    color = py_to_r(get_col(adata(), i))
                    title = isolate(input$color)
                    showlegend = FALSE
                    symbol = ~labels
                    symbols = all_symbols
                } else {
                    showNotification("Gene not found")
                    color = labels
                }
            }
            resubset(1) # split each subset into certain & uncertain subset or get back to original

            p <- plot_ly(
                x = x_emb_2d[, 1], y = x_emb_2d[, 2],
                text = text,
                color = color,
                #symbol = symbol,
                #symbols = symbols,  ## 30 shapes
                key = as.character(1:length(labels)),
                marker = list(size = isolate(input$dot_size)),
                type = 'scatter',
                mode = 'markers'
            )

            p <- p %>% config(modeBarButtonsToAdd = list(plot_options_btn))

            p <- p %>% layout(
                dragmode = "lasso",
                showlegend = showlegend,
                title = title,
                margin = list(t = 50))

            p <- theme_plot(p, theme_mode=isolate(input$theme_mode))

            p <- plotly_build(p)

            if (no_app == 0) {
                isolate(plotHistory(c(isolate(plotHistory()), list(p))))
                curPlot(length(plotHistory()))
            } else {
                tmp <- isolate(plotHistory())
                curPlot(length(isolate(plotHistory())))
                tmp[isolate(curPlot())] = list(p)
                isolate(plotHistory(tmp))
                navigate(1)
            }
        })
    })

    observeEvent(input$dot_size, {
        if (is_active(adata()) == FALSE) return()
        runjs(js.reset_marker_size)
    })

    observeEvent(input$plot_height, {
        if (is_active(adata()) == FALSE) return()
        runjs(js.reset_plot_height)
    })

    observeEvent(input$show_names, {
        if (is_active(adata()) == FALSE) return()
        no_append(1)
        replot(replot() + 1)
    })

    observeEvent(input$theme_mode, {
        if (is_active(adata()) == FALSE) return()
        runjs(js.reset_theme)
    })

    observeEvent(input$color, {
        if (is_active(adata()) == FALSE) return()
        replot(replot() + 1)
    })

    ###########################################################################
    # Plot Navigation

    toListenReplot <- reactive({
        list(curPlot(), navigate())
    })

    observeEvent(toListenReplot(), {
        if (curPlot() < 1) return()

        output$plot <- renderPlotly({
            isolate(navigate(0))
            p <- isolate(plotHistory())[[curPlot()]]
            p$x$layout$height = isolate(input$plot_height)
            for (i in seq_along(p$x$data))
               p$x$data[[i]]$marker$size = isolate(input$dot_size)
            p <- theme_plot(p, theme_mode=isolate(input$theme_mode))
            p <- p %>% toWebGL()
            return(p)
        })

        output$cell_names_outp <- renderUI({
            return(isolate(cellNamesTb())[[isolate(curPlot())]])
        })

        output$clustering_info <- renderUI({
            return(isolate(infoTb())[[isolate(curPlot())]])
        })
    })

    observeEvent(input$prevplot, {
        if (curPlot() < 1) return()

        if (curPlot() == 1) {
            showNotification("No more plots to show")
            return()
        }
        navigate(navigate() + 1)
        curPlot(curPlot() - 1)
    })

    observeEvent(input$nextplot, {
        if (curPlot() < 1) return()

        if (curPlot() == length(plotHistory())) {
            showNotification("This is the last plot")
            return()
        }
        curPlot(curPlot() + 1)
        navigate(navigate() + 1)
    })

    observeEvent(input$firstplot, {
        if (curPlot() < 1) return()
        curPlot(1)
        navigate(navigate() + 1)
    })

    observeEvent(input$lastplot, {
        if (curPlot() < 1) return()
        curPlot(length(plotHistory()))
        navigate(navigate() + 1)
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
