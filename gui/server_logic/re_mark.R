re_mark <- function(input, output, session, remark, pipe, deGenes, deButtons,
                    plotObj, rebutton, replot) {
    observe({ if (remark() > 0) {
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

        updateTabsetPanel(session, "switcher", selected = "DE")

        output$genes <- renderTable({
            # TODO change '0' to whatever the first element is
            # not all clustering algorithms return cluster ids starting with 0
            table <- data.frame(
                        pipe()$markers[['0']][['outp_names']], # gene_names
                        pipe()$markers[['0']][['diffs']], # logFC
                        pipe()$markers[['0']][['pvals']], # pvals
                        check.names = TRUE, # remove duplicates etc
                        stringsAsFactors = FALSE)
            table <- table[order(table[, 2], decreasing = TRUE),]
            colnames(table) <- c("DE genes", "logFC", "adj.pval")
            table[, 3] <- format(table[, 3], scientific = T)

            deGenes(table[, 1])
            return(table)
        },caption = tabletitle,
        caption.placement = getOption("xtable.caption.placement", "top"),
        caption.width = getOption("xtable.caption.width", NULL))

        output$DEbuttons <- renderUI({
            lapply(X = 1:length(deGenes()), FUN = function(i) {
                    actionButton(ns(paste(deGenes()[i], "_", seq = "")),
                                 deGenes()[i])})
        })

        rebutton(1)
        remark(0)
    }})

    observe({
    if (rebutton() > 0) {
        lapply(X = 1:length(deGenes()), FUN = function(i) {
            o <- observeEvent(input[[paste(deGenes()[i], "_", seq = "")]], {
                updateSelectInput(
                    session,
                    'color',
                    selected = deGenes()[i]
                )
                replot(1)
            })
            deButtons(c(deButtons(), isolate(o)))
        })
        rebutton(0)
    }})
}

reset_analysis_tabs <- function(output) {
    # Reset analysis tabs and switch to DE
    # By default the tabs will be erased only when we switch to them
    # so we first disable suspendWhenHidden, set them to null, and enable
    # it again
    outputOptions(output, "genes", suspendWhenHidden = FALSE)
    outputOptions(output, "DEbuttons", suspendWhenHidden = FALSE)
    outputOptions(output, "KEGG", suspendWhenHidden = FALSE)
    outputOptions(output, "GeneOntology", suspendWhenHidden = FALSE)
    outputOptions(output, "Markers", suspendWhenHidden = FALSE)
    outputOptions(output, "Msigdb", suspendWhenHidden = FALSE)
    output$DEbuttons = NULL
    output$genes = NULL
    output$GeneOntology = NULL
    output$KEGG = NULL
    output$Markers = NULL
    output$Msigdb = NULL
    outputOptions(output, "genes", suspendWhenHidden = TRUE)
    outputOptions(output, "DEbuttons", suspendWhenHidden = TRUE)
    outputOptions(output, "KEGG", suspendWhenHidden = TRUE)
    outputOptions(output, "GeneOntology", suspendWhenHidden = TRUE)
    outputOptions(output, "Markers", suspendWhenHidden = TRUE)
    outputOptions(output, "Msigdb", suspendWhenHidden = TRUE)
}
