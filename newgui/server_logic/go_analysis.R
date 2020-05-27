load('gui/obj/Hs.c2')

go_analysis <- function(input, output, session) {
    go_categories <- names(Hs.c5)
    len = length(go_categories)
    go_n <- integer(length = len)
    go_nde <- integer(length = len)
    go_pvals <- double(length = len)
    go_genes <- character(length = len)
    go_dispdat <- data.frame(go_categories, go_n, go_nde, go_pvals, go_genes,
                            stringsAsFactors = FALSE)

    colnames(go_dispdat) <- c("NAME","N","Intersect Length","P Value","Intersect Genes")
      output$GeneOntology <- renderTable({
        rownames(gene_ids_all)<-gene_ids_all[,1]
        withProgress(message = 'calculating Gene Ontology',detail=NULL, value = 0, {
          incProgress(1/3, detail = paste("Step: Getting gene IDs"))
          degenes_ids<-degenes_table_ord[1:input$mark_markers_n,1]
          degenes_int<-intersect(degenes_ids,rownames(gene_ids_all))
          degenes<-gene_ids_all[degenes_int,3]
          incProgress(2/3, detail = paste("Step: Calculating "))
          for (i in 1:nrow(go_dispdat)) {
            go_dispdat[i,2]<-length(Hs.c5[[i]])
            go_dispdat[i,3]<-length(intersect(degenes,Hs.c5[[i]]))
            go_dispdat[i,4]<-phyper(length(intersect(degenes,Hs.c5[[i]])),length(Hs.c5[[i]]),ncol(scdata_subset)-1-length(Hs.c5[[i]]),length(degenes),lower.tail = F)
            int_genes_ids<-intersect(degenes,Hs.c5[[i]])
            rownames(gene_ids_all)<-gene_ids_all[,3]
            int_genes<-gene_ids_all[int_genes_ids,1]
            if (length(int_genes)>0){
              go_dispdat[i,5]<-(paste(int_genes,collapse=", "))
            } else {
              go_dispdat[i,5]<-as.character("0")
            }
          }
          incProgress(3/3, detail = paste("Step: Getting Geneontology"))
          go_filt<-go_dispdat[which(go_dispdat[,4]<0.05),]
          go_ord<-go_filt[order(go_filt[,4]),]
          showNotification("GeneOntology calculation finished")
          output$downloadGO <- downloadHandler(
            filename = function() {
              paste("GO_data", ".csv", sep = "")
            },
            content = function(file) {
              write.csv(go_ord, file, row.names = FALSE)
            }
          )
          go_ord[,4]<-format(go_ord[,4],scientific=T)
          head(go_ord,n = 10)
        })
      },bordered = T)
}

