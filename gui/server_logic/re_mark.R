re_mark <- function(input, output, session, remark, pipe, deGenes, deButtons,
                    rebutton, replot) {
    observe({ if (remark() > 0) {
        rebutton(rebutton() + 1)
        remark(0)
        updateTabsetPanel(session, "switcher", selected = "DE")
        reset_analysis_tabs(output)
        if (length(deButtons()) > 0)
            for (i in 1:length(deButtons()))
                deButtons()[[i]]$destroy()
        deButtons(c())
        deGenes(c())

        ns <- session$ns

        #reset_analysis_tabs(output)
        s1=as.character(input$subset1)
        s2=as.character(input$subset2)
        tabletitle=paste(s1," (vs. ",s2,")",sep="")
        output$titleDE <- renderText(tabletitle)
        updateTabsetPanel(session, "switcher", selected = "DE")

        shinyInput <- function(FUN, len, id, ns, ...) {  ## for inputting buttons
            inputs <- character(len)
            for (i in seq_len(len)) {
                inputs[i] <- as.character(FUN(paste0(id, i), ...))
            }
            inputs
        }

        output$DEtable <- DT::renderDataTable({
            # TODO change '0' to whatever the first element is
            # not all clustering algorithms return cluster ids starting with 0
            table <- data.frame(
                        pipe()$markers[['0']][['outp_names']], # gene_names
                        pipe()$markers[['0']][['diffs']], # logFC
                        pipe()$markers[['0']][['pvals']], # pvals
                        Actions = shinyInput(actionButton, length(pipe()$markers[['0']][['outp_names']]),
                                             'button_',
                                             label = "Show Expression",
                                             #onclick = paste0('Shiny.onInputChange(\"' , ns("select_button"), '\", this.id)')
                                             onclick = sprintf("Shiny.onInputChange('%s', this.id)", ns("select_button")
                                             ),
                        ),
                        check.names = TRUE, # remove duplicates etc
                        stringsAsFactors = FALSE)
            table <- table[order(table[, 2], decreasing = TRUE),]
            colnames(table) <- c("DE genes", "logFC", "adj.pval","Show Expression Level")
            table[, 3] <- format(table[, 3], scientific = T)

            deGenes(table[, 1])
            output$downloadDE <- downloadHandler(
                filename = function() {
                    paste0(pipe()$dataset, "_DE_genes", ".csv")
                },
                content = function(file) {
                    write.csv(table, file, row.names = FALSE)
                }
            )

            return(table)
        },
        #caption = tabletitle,
        #caption.placement = getOption("xtable.caption.placement", "top"),
        #caption.width = getOption("xtable.caption.width", NULL)
        server = FALSE
        , escape = FALSE
        ,selection = 'none'
        )

        o<-observeEvent(input$select_button, {
            if (identical(input$select_button,'reset'))
            {
                return()
            }
            selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
            #showNotification(as.character(selectedRow))
            if (length(pipe()$markers[['0']][['outp_names']])>=selectedRow)
            {
                selected_gene=pipe()$markers[['0']][['outp_names']][[selectedRow]]
                updateSelectInput(
                    session,
                    'color',
                    selected = selected_gene
                )
            }
            #myValue$employee <<- paste('click on ',df$data[selectedRow,1])
        })
        deButtons(c(deButtons(),o))
        # output$DEbuttons <- renderUI({
        #     lapply(X = 1:length(deGenes()), FUN = function(i) {
        #             actionButton(ns(paste(deGenes()[i], "_", seq = "")),
        #                          deGenes()[i])})
        # })
    }})

    # observe({
    # if (rebutton() > 0) {
    #     lapply(X = 1:length(deGenes()), FUN = function(i) {
    #         o <- observeEvent(input[[paste(deGenes()[i], "_", seq = "")]], {
    #             updateSelectInput(
    #                 session,
    #                 'color',
    #                 selected = deGenes()[i]
    #             )
    #         })
    #         deButtons(c(deButtons(), isolate(o)))
    #     })
    #     rebutton(0)
    # }})
}

reset_analysis_tabs <- function(output) {
    # Reset analysis tabs and switch to DE
    # By default the tabs will be erased only when we switch to them
    # so we first disable suspendWhenHidden, set them to null, and enable
    # it again
    outputOptions(output, "DEtable", suspendWhenHidden = FALSE)
    outputOptions(output, "DEbuttons", suspendWhenHidden = FALSE)
    outputOptions(output, "KEGGtable", suspendWhenHidden = FALSE)
    outputOptions(output, "GOtable", suspendWhenHidden = FALSE)
    outputOptions(output, "Markerstable", suspendWhenHidden = FALSE)
    outputOptions(output, "MSIGDBtable", suspendWhenHidden = FALSE)
    outputOptions(output, "Markerstableuser", suspendWhenHidden = FALSE)
    output$DEbuttons = NULL
    output$DEtable = NULL
    output$GOtable = NULL
    output$KEGGtable = NULL
    output$Markerstable = NULL
    output$Markerstableuser = NULL
    output$MSIGDBtable = NULL
    outputOptions(output, "DEtable", suspendWhenHidden = TRUE)
    outputOptions(output, "DEbuttons", suspendWhenHidden = TRUE)
    outputOptions(output, "KEGGtable", suspendWhenHidden = TRUE)
    outputOptions(output, "GOtable", suspendWhenHidden = TRUE)
    outputOptions(output, "Markerstable", suspendWhenHidden = TRUE)
    outputOptions(output, "MSIGDBtable", suspendWhenHidden = TRUE)
    outputOptions(output, "Markerstableuser", suspendWhenHidden = TRUE)
}
