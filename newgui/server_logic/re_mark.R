re_mark <- function(input, output, session, remark, pipe, deGenes, deButtons) {
    observe({ if (remark() > 0) {
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
            deGenes(table[1:length(pipe()$markers[['0']][['outp_names']]), 1])
            return(table)
        })

        if (length(deButtons()) > 0)
            for (i in 1:length(deButtons()))
                deButtons()[[i]]$destroy()
        deButtons(c())

        output$DEbuttons <- renderUI({
            lapply(X = 1:length(deGenes()), FUN = function(i) {
                    actionButton(paste(deGenes()[i], " ", seq = ""),
                                 paste(deGenes()[i], " ", seq = ""))})
        })


        lapply(X = 1:length(deGenes()), FUN = function(i) {
            o <- observeEvent(input[[paste(deGenes()[i], " ", seq = "")]], {
                showNotification(
                    paste("Showing ", deGenes()[i], "'s expression", sep = ""),
                    duration=5)
            })
            deButtons(c(deButtons(), isolate(o)))
        })
        remark(0)
    }})
}
