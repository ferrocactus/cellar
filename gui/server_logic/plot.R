source("gui/server_logic/plot_options.R")

# Define the number of colors you want
plot <- function(input, output, session, replot, adata, activeDataset,
                 setNames, setPts, reset, resubset, reinfo, relabel,
                 retheme, info_val) {
    ns <- session$ns
    plot_index <- reactiveVal(0)
    plot_count <- reactiveVal(0)
    #lider_update <- reactiveVal(0)
    main_plot_val <- reactiveVal(NULL)

    # triggers when replot is set to 1
    observeEvent(replot(), {
        if (replot() < 1) return()
        isolate(replot(0))

        req(adata())
        if (!py_has_attr(adata()$obs, 'labels')) return()

        relabel(relabel() + 1)
        reinfo(reinfo() + 1)

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
            if (isolate(input$show_names) == 'show_names' && isolate(input$color) == 'Clusters')
                color = paste0(labels, ": ", label_names)
            else if (isolate(input$color) == 'Clusters')
                color = labels
            else if (isolate(input$color) == 'Uncertainty') {
                if (anyNA(as.integer(isolate(input$n_neighbors))) == TRUE) {
                    n_neighbors = as.integer(sqrt(length(labels)))
                } else if (as.integer(isolate(input$n_neighbors)) > as.integer(length(labels) / 2)) {
                    n_neighbors = as.integer(sqrt(length(labels)))
                } else {
                    n_neighbors = as.integer(isolate(input$n_neighbors))
                }

                showNotification("Calculating Uncertainty")
                withProgress(message = "Please Wait", value = 0, {
                    incProgress(1 / 3, detail = "Data processing...")
                    cellar$safe(cellar$get_neighbors,
                                x = adata(),
                                n_neighbors = n_neighbors)
                    incProgress(1 / 3, detail = "Calculating uncertainty...")
                    cellar$safe(cellar$uncertainty,
                                x = adata(),
                                n_neighbors = n_neighbors)
                })
                showNotification("Finished")

                color = as.numeric(py_to_r(adata()$obs['uncertainty']))
                title = "Uncertainty"
                showlegend = FALSE
                symbol = ~labels
                symbols = all_symbols
                cell_ids = as.character(py_to_r(adata()$obs$index$to_numpy()))
                certainty = py_to_r(adata()$uns['uncertainty_text'])
                text = ~paste0(replicate(length(label_names), "Label: "),
                            label_names, replicate(length(label_names), ", Uncertainty: "),
                            as.character(py_to_r(adata()$obs['uncertainty'])), '\n', as.character(certainty))
            } else {
                # show gene expression level:
                gene_names = py_to_r(get_all_gene_names(adata()))
                i = which(gene_names == isolate(input$color))[1]
                if (!is.null(i)) {
                    
                    text = ~paste("Label: ", label_names,'\nColor value:',as.character((color))) # text shows expression value
                    color = py_to_r(get_col(adata(), i))
                    
                    color=color-min(color)+1 # make min=1
                    color = log(color)  # min=0
                    
                     
                    if (identical(NULL,input$value_t)==FALSE){
                        t=as.numeric(input$value_t)  # threshold
                        if (t>0){
                            for (i in 1:length(color)){
                                if (color[i]<t){
                                    color[i]=0  # 'gray'
                                }
                            }
                        }
                    }
                   
                    # if (is.na(min(color[color > 0]))==FALSE){
                    #     mini=min(color[color > 0])
                    #     for (i in 1:length(color)){
                    #         if (color[i]==0){
                    #             color[i]=mini  #
                    #         }
                    #     }
                    # }
                    
                    
                    
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
                mode = 'markers',
                height = isolate(input$plot_height),
                source = 'M'
            )


            p <- p %>% config(modeBarButtonsToAdd = list(plot_options_btn))

            p <- p %>% layout(
                dragmode = "lasso",
                showlegend = showlegend,
                title = title,
                margin = list(t = 50))

            p <- theme_plot(p, theme_mode = isolate(input$theme_mode))
            p <- plotly_build(p)

            isolate(main_plot_val(p))
            p <- p %>% toWebGL()

            output$plot <- renderPlotly({ p })
        })
    })
    
    observeEvent(input$color,{
        if (input$color!='Uncertainty' && input$color != 'Clusters'){
            gene_names = py_to_r(get_all_gene_names(adata()))
            i = which(gene_names == isolate(input$color))[1]
            color = py_to_r(get_col(adata(), i))
            color=color-min(color)+1 # make min=1
            color = log(color)  # min=0
           
            output$threshold_slider <- renderUI({
                sliderInput(
                    ns("value_t"),
                    label="Value threshold",
                    min=0,max=signif(max(color),digits=3),
                    step=(max(color)/30),
                    value=0
                )
            })
        }
        else{
            output$threshold_slider <- renderUI({
               
            })
        }
    })
    
    observe({
        req(info_val$cellNames)
        output$cell_names_outp <- renderUI({ info_val$cellNames })
        req(info_val$configs)
        output$clustering_info <- renderUI({ info_val$configs })
    })

    observeEvent(input$value_t, {
        req(adata())
        if (is.na(as.numeric(input$value_t))==FALSE){
            if ((input$color != 'Clusters') && (input$color != 'Uncertainty')){
                replot(replot() + 1)
            }
        }
    })
    
    observeEvent(input$dot_size, {
        req(adata())
        runjs(js.reset_marker_size)
    })

    observeEvent(input$plot_height, {
        req(adata())
        runjs(js.reset_plot_height)
    })

    observeEvent(input$show_names, {
        req(adata())
        replot(replot() + 1)
    })

    observeEvent(input$theme_mode, {
        req(adata())
        runjs(js.reset_theme)
    })

    observeEvent(input$color, {
        req(adata())
        replot(replot() + 1)
    })

    observeEvent(input$store_plot, {
        req(main_plot_val())

        if (plot_count() == 5) {
            showNotification("Cannot have more than 5 stored plots.")
            return()
        }

        isolate(plot_index(isolate(plot_index()) + 1))
        isolate(plot_count(isolate(plot_count()) + 1))

        plot_i = as.character(isolate(plot_index()))
        title = paste0("Plot ", plot_i, " (", activeDataset(), ")")
        plot_id = paste0("plot", plot_i)
        collapse_btn_id = paste0("collapse_cell_names", plot_i)
        configs_id = paste0("clustering_info", plot_i)
        cell_names_id = paste0("cell_names_outp", plot_i)
        observer_id = paste0("observer", plot_i)

        appendTab(
            "tabset", tabPanel(
                title,
                plotlyOutput(ns(plot_id), height = "100%"),
                div(
                    class = "cell_names_div",
                    list(
                        actionButton(
                            ns(collapse_btn_id),
                            "View additional info",
                            class = "collapsebtn"),
                        splitLayout(
                            htmlOutput(ns(configs_id)),
                            htmlOutput(ns(cell_names_id)))))),
            select = TRUE)

        observeEvent(input[[collapse_btn_id]], {
            shinyjs::toggle(cell_names_id)
            shinyjs::toggle(configs_id)
        })

        output[[cell_names_id]] <- renderUI({ isolate(info_val$cellNames) })

        output[[configs_id]] <- renderUI({ isolate(info_val$configs) })

        output[[plot_id]] <- renderPlotly({
            p <- isolate(main_plot_val())
            p$x$layout$height = isolate(input$plot_height)
            for (i in seq_along(p$x$data))
                p$x$data[[i]]$marker$size = isolate(input$dot_size)
            p <- theme_plot(p, theme_mode = isolate(input$theme_mode))
            p <- p %>% toWebGL()
            p$x$layout$title <- title
            p$x$source <- 'P'
            retheme(1)
            return(p)
        })
    })

    observeEvent(input$delete_plot, {
        if (isolate(input$tabset) == "Main Plot") {
            showNotification("Cannote delete main plot.")
            return()
        }

        isolate(plot_count(isolate(plot_count()) - 1))
        removeTab("tabset", isolate(input$tabset))
    })
    ###########################################################################

    # Store selected cells
    observeEvent(input$store_lasso, {
        req(main_plot_val())
        if (!py_has_attr(adata()$obs, 'labels')) return()
        req(input$newsubset)

        if (isolate(input$tabset) != "Main Plot") {
            showNotification("You can only select in the main plot.")
            return()
        }

        if (as.character(isolate(input$newsubset)) %in% setNames()) {
            showNotification("Name already exists")
            return()
        }

        if (substr(as.character(isolate(input$newsubset)), 1, 7) == "Cluster") {
            showNotification("Reserved name, please choose another.")
            return()
        }

        d <- event_data("plotly_selected", source='M')
        req(d$key)

        keys <- as.numeric(d$key)
        msg <- cellar$safe(store_subset,
            adata = adata(),
            name = isolate(input$newsubset),
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
            extension <- tolower(isolate(input$plot_download_format))
            paste0(tools::file_path_sans_ext(basename(activeDataset())),
                   "_plot.", extension)
        },
        content = function(fname) {
            req(main_plot_val())
            extension <- tolower(isolate(input$plot_download_format))
            if (extension == 'html') {
                htmlwidgets::saveWidget(
                    as_widget(main_plot_val()), fname,
                    selfcontained = TRUE)
            } else {
                withProgress(message = "Rendering plot", value = 0, {
                    incProgress(1 / 2)
                    withr::with_dir(dirname(fname),
                                    orca(main_plot_val(), basename(fname),
                                    format = extension))
                })
            }
        }
    )
    return(main_plot_val)
}
