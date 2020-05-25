      observeEvent(input$ssclurun, {
        if(exists("updated_new_labels")) {
          tryCatch({
            df <- runSSClu(pipe, updated_new_labels, input)
          }, error = function(e) {
            df <- "An error occurred."
          })

          if (is.character(df) & length(df) == 1) {
            showNotification(df)
          } else {
            markers <- pipe$markers


            #get the gene expression data
            expr_data = data.frame(
              matrix(pipe$x, ncol = length(pipe$col_ids),
                     dimnames = list(1:length(pipe$x[,1]), pipe$col_ids))
            )
            expr_data$cluster = df$y
            ##read in the marker from JSON
            hypergeom <- getHypergeom("markers/cell_type_marker.json")
            markers_genelists_list <- getMarkerGeneList(
              "markers/cell_type_marker.json")


            newlabs <- df[, 3]
            names(newlabs) <- rownames(df)
            labeldats <- levels(as.factor(df[,3]))

            output$plot <- renderPlotly({
              #factorize cluster labels (discrete instead of continuous)
              if (input$color == "cluster") {
                plotcols = as.factor(expr_data[[input$color]])
              } else {
                plotcols = expr_data[[input$color]]
              }
              plot_ly(
                df, x = df$x1, y = df$x2,
                text = ~paste("label: ", as.factor(df$y)),
                color = plotcols,
                key = row.names(df),
                type = 'scatter',
                mode = 'markers'
              ) %>% layout(dragmode = "lasso",
                           title = paste("Value of ", input$color, sep=""))
            })

            ### updated plot
            observeEvent(input$labelupd, {
              #d <- event_data("plotly_selected")
              keys<-s1
              keysNotification(as.character(input$newlabels))

              newlabs[keys]<<-as.character(input$newlabels)


              assign("updated_new_labels", newlabs, envir = env)
              output$Plot2 <- renderPlotly({
                plot_ly(
                  df, x = df[, 1], y = df[, 2],
                  text = ~paste("label: ", as.factor(newlabs)),
                  color = as.factor(newlabs),
                  type = 'scatter',
                  mode = 'markers'
                ) %>% layout(dragmode = "lasso",
                             title = paste("Value of ", input$labelupd, sep=""))
              })
              output$downlabels <- downloadHandler(
                filename = function() {
                  paste("Updated_labels", ".csv", sep = "")
                },
                content = function(file) {
                  write.csv(newlabs, file, row.names = FALSE)
                }
              )
            })
          }
        }
      })
