source("gui/server_logic/plot_options.R")

EPS = 0.01

# Define the number of colors you want
plot <- function(input, output, session, replot, adata, activeDataset,
                 setNames, setPts, reset, resubset, reinfo, relabel,
                 retheme, info_val_cellNames, info_val_configs) {
  ns <- session$ns
  plot_index <- reactiveVal(0)
  plot_count <- reactiveVal(0)
  double_plot <- reactiveVal(FALSE)
  #double_plot2 <- reactiveVal(FALSE)
  plot_cell_labels <- list()
  plot_cell_labels_split <- NULL
  #plot_cell_labels_split2 <- NULL
  GRAY <- c(220, 220, 220)
  #lider_update <- reactiveVal(0)
  main_plot_val <- reactiveVal(NULL)
  #main_plot_val2 <- reactiveVal(NULL)
  output$threshold_slider <- NULL
  outputOptions(output, "threshold_slider", suspendWhenHidden = FALSE)

  trigger_threshold <- reactiveVal(FALSE)

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
      colors = NULL

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
        print(n_neighbors)

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
      } else if (isolate(input$color2) == "None") {
        # show gene expression level:
        gene_names = py_to_r(get_all_gene_names(adata()))
        i = which(gene_names == isolate(input$color))[1]
        if (!is.null(i)) {
          color = py_to_r(get_col(adata(), i))
          text = ~paste("Label: ", label_names,'\nColor value:',as.character((color)))

          m = signif(min(color) - EPS, digits=3)
          M = signif(max(color) + EPS, digits=3)

          if (isolate(trigger_threshold()) == TRUE) {
            v1 = isolate(input$value_t)[1]
            v2 = isolate(input$value_t)[2]
          } else {
            v1 = m
            v2 = M
          }

          colors <- function(vals) {
            c_func = colorRamp(c("#440154", "#27828D", "#FDE725"))

            if (!is.null(v1)) { # New gene
              # Need to be in [0, 1] range so divide by max
              min_t = as.numeric(v1) / M
              max_t = as.numeric(v2) / M
              color_matrix = matrix(0, length(vals), 3)
              for (i in 1:length(vals)) {
                if (vals[i] < min_t) color_matrix[i, 1:3] = GRAY
                else if (vals[i] >= max_t) color_matrix[i, 1:3] = GRAY
                # Otherwise map [min_t, max_t] to [0, 1]
                else color_matrix[i, 1:3] = c_func((vals[i] - min_t) / (max_t - min_t))
              }
              return(color_matrix)
            } else { # default
              return(c_func(vals))
            }
          }
          title = isolate(input$color)
          showlegend = FALSE
        } else {
          showNotification("Gene not found")
          color = labels
        }
      } else {
        # Co-expression level
        gene_names = py_to_r(get_all_gene_names(adata()))
        gene1 = which(gene_names == isolate(input$color))[1]
        gene2 = which(gene_names == isolate(input$color2))[1]
        color1 = py_to_r(get_col(adata(), gene1)) # true expression
        color2 = py_to_r(get_col(adata(), gene2)) # true expression
        text = ~paste0("Label: ", label_names,
                       '\n', input$color, ': ', as.character((signif(color1, digits=3))),
                       '\n', input$color2, ': ', as.character((signif(color2, digits=3))))

        m1 = signif(min(color1) - EPS, digits=3)
        M1 = signif(max(color1) + EPS, digits=3)
        m2 = signif(min(color2) - EPS, digits=3)
        M2 = signif(max(color2) + EPS, digits=3)

        # Map color1 to [0, 1]
        color1 = (color1 - m1) / (M1 - m1)
        # Map color2 to [0, 1]
        color2 = (color2 - m2) / (M2 - m2)

        if (isolate(trigger_threshold()) == TRUE) {
          # slider bar thresholds mapped to [0, 1]
          v11 = (isolate(input$value_t)[1] - m1) / (M1 - m1)
          v12 = (isolate(input$value_t)[2] - m1) / (M1 - m1)
          v21 = (isolate(input$value_t_2)[1] - m2) / (M2 - m2)
          v22 = (isolate(input$value_t_2)[2] - m2) / (M2 - m2)
        } else {
          v11 = 0 - EPS
          v12 = 1 + EPS
          v21 = 0 - EPS
          v22 = 1 + EPS
        }

        for (i in 1:length(color1)) {
          if (color1[i] < v11 || color1[i] > v12) color1[i] = 0
          # color1 should have the same length as color2
          if (color2[i] < v21 || color2[i] > v22) color2[i] = 0
        }

        color = pmin(color1, color2)

        colors <- function(vals) {
          c_func = colorRamp(c("#440154", "#27828D", "#FDE725"))

          color_matrix = matrix(0, length(vals), 3)
          for (i in 1:length(vals)) {
            if (vals[i] < EPS) color_matrix[i, 1:3] = GRAY
            else {
              # all non-outlier vals are in [0.5, 1]
              #v = (vals[i] - 0.5) * 2
              #v = max(min(1, v), 0) # ensure v is in [0, 1]
              color_matrix[i, 1:3] = c_func(vals[i])
            }
          }
          return(color_matrix)
        }

        title = paste(isolate(input$color), 'co.', isolate(input$color2))
        showlegend = FALSE
      }
      resubset(1) # split each subset into certain & uncertain subset or get back to original

      p <- plot_ly(
        x = x_emb_2d[, 1], y = x_emb_2d[, 2],
        text = text,
        color = color,
        colors = colors,
        #symbol = symbol,
        #symbols = symbols,  ## 30 shapes
        key = as.character(1:length(labels)),
        marker = list(size = isolate(input$dot_size)),
        type = 'scatter',
        mode = 'markers',
        height = isolate(input$plot_height),
        source = 'M'
      )

      p2 <- plot_ly(
        x = x_emb_2d[, 1], y = x_emb_2d[, 2],
        text = text,
        color = color,
        colors = colors,
        #symbol = symbol,
        #symbols = symbols,  ## 30 shapes
        key = as.character(1:length(labels)),
        marker = list(size = isolate(input$dot_size)),
        type = 'scatter',
        mode = 'markers',
        #height = isolate(input$plot_height),
        source = 'M'
      )

      p <- p %>% config(
        displaylogo = FALSE,
        #modeBarButtonsToAdd = list(plot_options_btn),
        displayModeBar = TRUE)

      # p2 <- p2 %>% config(
      #   displaylogo = FALSE,
      #   #modeBarButtonsToAdd = list(plot_options_btn),
      #   displayModeBar = TRUE)

      p <- p %>% layout(
        dragmode = "lasso",
        showlegend = showlegend,
        title = title,
        margin = list(t = 50))

      p <- theme_plot(p, theme_mode = isolate(input$theme_mode))
      p <- plotly_build(p)

      isolate(main_plot_val(p))
      p <- p %>% toWebGL()

      # p2 <- p2 %>% layout(
      #   dragmode = "lasso",
      #   showlegend = showlegend,
      #   title = title,
      #   margin = list(t = 50))

      p2 <- theme_plot(p2, theme_mode = isolate(input$theme_mode))
      #p2 <- plotly_build(p2)

      isolate(main_plot_val(p))
      #isolate(main_plot_val2(p2))
      #p2 <- p2 %>% toWebGL()

      output$plot <- renderPlotly({ p })
      #output$plotg <- renderPlotly({ p2 })#,height=input$plot_height)
      #output$plotg2 <- NULL#
    })
  })

  observeEvent(input$color, {
    if (input$color!='Uncertainty' && input$color != 'Clusters'){
      gene_names = py_to_r(get_all_gene_names(adata()))
      i = which(gene_names == isolate(input$color))[1]
      color = py_to_r(get_col(adata(), i))

      m = signif(min(color) - EPS, digits=3)
      M = signif(max(color) + EPS, digits=3)
      step = signif((M - m) / 30, digits=3)

      isolate(trigger_threshold(FALSE))
      updateSliderInput(
        session = session,
        inputId = "value_t",
        min = m, max = M,
        value = c(m, M),
        step = step
      )
    } else {
      isolate(trigger_threshold(FALSE))
      updateSliderInput(
        session = session,
        inputId = "value_t",
        min = 0, max = 1,
        value = c(0, 1),
        step = 1
      )
    }
  })

  observeEvent(input$color2, {
    if (input$color!='Uncertainty' && input$color != 'Clusters' && input$color2 != 'None') {
      gene_names = py_to_r(get_all_gene_names(adata()))
      i = which(gene_names == isolate(input$color2))[1]
      color = py_to_r(get_col(adata(), i))

      m = signif(min(color) - EPS, digits=3)
      M = signif(max(color) + EPS, digits=3)
      step = signif((M - m) / 30, digits=3)

      if (input$color2 != 'None')
        isolate(trigger_threshold(FALSE))
      updateSliderInput(
        session = session,
        inputId = "value_t_2",
        min = m, max = M,
        value = c(m, M),
        step = step
      )
    } else {
      if (input$color2 != 'None')
        isolate(trigger_threshold(FALSE))
      updateSliderInput(
        session = session,
        inputId = "value_t_2",
        min = 0, max = 1,
        value = c(0, 1),
        step = 1
      )
    }
  })

    observeEvent(info_val_cellNames(), {
        output$cell_names_outp <- renderUI({
            req(info_val_cellNames())
            return(info_val_cellNames())
        })
    })

    observeEvent(info_val_configs(), {
        output$clustering_info <- renderUI({
            req(info_val_configs())
            return(info_val_configs())
        })
    })

  observeEvent(input$value_t, {
    if (trigger_threshold() == FALSE) {
      isolate(trigger_threshold(TRUE))
      return()
    }
    req(adata())
    if (!is.null(input$value_t[0])){
      if ((input$color != 'Clusters') && (input$color != 'Uncertainty')){
        replot(replot() + 1)
      }
    }
  })

  observeEvent(input$value_t_2, {
    if (trigger_threshold() == FALSE) {
      isolate(trigger_threshold(TRUE))
      return()
    }
    if (isolate(input$color2) == 'None') return()
    req(adata())
    if (!is.null(input$value_t[0])){
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

  observeEvent(input$color2, {
    req(adata())
    if (isolate(input$color) == 'Clusters' || isolate(input$color) == 'Uncertainty') return()
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
    transfer_labels_id = paste0("transfer_labels", plot_i)

    # Store cell ids and labels
    cell_names = py_to_r(get_obs_names(adata()))
    labels = py_to_r(get_label_names(adata()))
    plot_cell_labels[[transfer_labels_id]] = list(cell_names, labels)

    appendTab(
      "tabset", tabPanel(
        title,
        plotlyOutput(ns(plot_id), height = "100%"),
        div(
          class = "cell_names_div",
          list(
            actionButton(
              ns(transfer_labels_id),
              "Move labels to Main Plot"),
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

    # observeEvent(input[['collapse_cell_names2']], {
    #   shinyjs::toggle('clustering_info2')
    #   shinyjs::toggle('cell_names_outp2')
    # })

    observeEvent(input[[transfer_labels_id]], {
      lbs = plot_cell_labels[[transfer_labels_id]]

      match_labels(adata(), lbs[[1]], lbs[[2]])

      replot(replot() + 1)
      reset(reset() + 1)
      resubset(resubset() + 1)
    })

    output[[cell_names_id]] <- renderUI({ isolate(info_val_cellNames()) })

    output[[configs_id]] <- renderUI({ isolate(info_val_configs()) })

    output[['cell_names_outp2']] <- renderUI({ isolate(info_val_cellNames()) })

    output[['clustering_info2']] <- renderUI({ isolate(info_val_configs()) })

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

  insertUI(
    selector="#ns-plots",
    where="afterEnd",
    ui=plotlyOutput(ns("plot"), height="100%")
  )

  observeEvent(input$delete_plot, {
    if (isolate(input$tabset) == "Main Plot") {
      showNotification("Cannote delete main plot.")
      return()
    }

    isolate(plot_count(isolate(plot_count()) - 1))
    removeTab("tabset", isolate(input$tabset))
  })

  # observeEvent(input$split_plot2, {
  #   req(main_plot_val2())
  #   if (double_plot2() == TRUE) {
  #     output$plotg2=NULL
  #     isolate(double_plot2(FALSE))
  #     plot_cell_labels_split2 <- NULL
  #   }
  #   else {
  #           output$plotg2 <- renderPlotly({
  #             p = isolate(main_plot_val2())
  #             return(p)
  #           })
  #
  #           isolate(double_plot2(TRUE))
  #         }
  #     })

  observeEvent(input$split_plot, {
    req(main_plot_val())

    if (double_plot() == TRUE) {
      removeUI(selector="#ns-plot")
      removeUI(selector="#ns-plot2")
      removeUI(selector="#ns-double_split")
      removeUI(selector="#ns-move_labels_split")
      insertUI(
        selector="#ns-plots",
        where="beforeBegin",
        ui=plotlyOutput(ns("plot"), height="100%")
      )
      output$plot <- renderPlotly({
        p = isolate(main_plot_val())
        p <- p %>% toWebGL()
        return(p)
      })
      plot_cell_labels_split <- NULL
      isolate(double_plot(FALSE))
    } else {
      # Store cell ids and labels
      cell_names = py_to_r(get_obs_names(adata()))
      labels = py_to_r(get_label_names(adata()))
      plot_cell_labels_split = list(cell_names, labels)

      removeUI(selector="#ns-plot")
      insertUI(
        selector="#ns-plots",
        where="beforeBegin",
        ui=splitLayout(
          id=ns("double_split"),
          plotlyOutput(ns("plot"), height="100%"),
          plotlyOutput(ns("plot2"), height="100%")
        )
      )
      insertUI(
        selector="#ns-split_plot",
        where="afterEnd",
        ui=actionButton(
          ns("move_labels_split"),
          "Move labels"
        )
      )

      observeEvent(input$move_labels_split, {
        lbs = plot_cell_labels_split

        match_labels(adata(), lbs[[1]], lbs[[2]])

        replot(replot() + 1)
        reset(reset() + 1)
        resubset(resubset() + 1)
      })

      output$plot <- renderPlotly({
        p = isolate(main_plot_val())
        p <- p %>% toWebGL()
        return(p)
      })

      output$plot2 <- renderPlotly({
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
      isolate(double_plot(TRUE))
    }
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

print_app <- function(widget) {
  temp <- paste(tempfile('plotly'), 'html', sep = '.')
  htmlwidgets::saveWidget(widget, temp, selfcontained = TRUE)

  system(paste0("firefox ", temp))
}
