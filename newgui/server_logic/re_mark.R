re_mark <- function(input, output, session, remark, pipe) {
    observe({ if (remark() > 0) {
        output$genes <- renderTable({
            table <- data.frame(
                        pipe()$markers[['0']]['outp_names'], # gene_names
                        pipe()$markers[['0']]['diffs'], # logFC
                        pipe()$markers[['0']]['pvals'], # pvals
                        check.names = TRUE, # remove duplicates etc
                        stringsAsFactors = FALSE)
            table <- table[order(table[, 2], decreasing = TRUE),]
            print('gor here')
            colnames(table) <- c("DE genes", "logFC", "adj.pval")
            return(table)
        })

        #DEgenes<-degenes_table_ord[1:input$mark_markers_n,1] # a vector of characters

        #for (i in 1:length(DEgenes)) {
        #    assign("degenenames",c(degenenames,paste(DEgenes[i]," ",seq="")),envir=env)
        #}
        #  output$DEbuttons <- renderUI({
        #    lapply(
        #      X = 1:length(DEgenes),
        #      FUN = function(i){
        #        actionButton(paste(DEgenes[i]," ",seq=""),paste(DEgenes[i]," ",seq=""))
        #      }
        #    )
        #  })



          #selecteddat<-scdata_subset[as.numeric(keys),1:ncol(scdata_subset)-1]
          #restdat<-scdata_subset[-as.numeric(keys),1:ncol(scdata_subset)-1]

    }})
}
