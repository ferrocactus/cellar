library(shiny)
library(stringr)

source("gui/sidebar/options.R")

server <- shinyServer(function(input, output, session) {
    observeEvent(input$dim_method, {
        if(input$dim_method == "PCA"){
            shinyjs::enable(selector = "[type=radio][value=pca_auto]")
        } else {
            shinyjs::disable(selector = "[type=radio][value=pca_auto]")
            updateRadioButtons(
                session,
                "dim_options",
                selected = "pca_manual"
            )
        }
    })

    observeEvent(input$dim_options, {
        if(input$dim_options == "pca_manual"){
            shinyjs::enable('dim_n_components')
        } else {
            shinyjs::disable('dim_n_components')
        }
    })

    observeEvent(input$clu_method, {
        if(input$clu_method %in% options$clu_no_n_clusters){
            shinyjs::disable('clu_n_clusters')
            shinyjs::disable('eval_method')
        } else {
            shinyjs::enable('clu_n_clusters')
            shinyjs::enable('eval_method')
        }
    })

    observeEvent(input$ensemble_checkbox, {
        if('all' %in% str_split(input$ensemble_checkbox, " ")){
            for (x in options$clu_ensemble) {
                if (x != 'all') {
                    shinyjs::disable(
                        selector = paste("#ensemble_checkbox input[value='",
                                    x, "']", sep="")
                    )
                }
            }
        } else {
            for (x in options$clu_ensemble) {
                shinyjs::enable(
                    selector = paste("#ensemble_checkbox input[value='",x, "']",
                                     sep="")
                )
            }
        }
    })
})
