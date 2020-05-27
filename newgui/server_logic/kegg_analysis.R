#load('gui/obj/kegg_genelists')
load('gui/obj/Go_kegg_lists')
load('gui/obj/keggidtoname')
load('gui/obj/gene_ids_all')

kegg_analysis <- function(input, output, session, deGenes, pipe) {
    kegg_categories <- kegg_id_toname[names(kegg_genelists)]
    len <- length(kegg_categories)
    kegg_dispdat <- data.frame(kegg_categories,
                               integer(length = len),
                               integer(length = len),
                               double(length = len),
                               character(length = len),
                               stringsAsFactors = FALSE)
    colnames(kegg_dispdat) <- c("Name", "n", "Intersection Length",
                                "pval", "Intersection Genes")

    observe({
    if (length(deGenes()) > 0) {
        nc = length(pipe()$col_ids)

        output$KEGG <- renderTable({
            rownames(gene_ids_all) <- gene_ids_all[, 1]

            withProgress(message = 'calculating KEGG', detail = NULL,
                         value = 0, {
                incProgress(1/3, detail = paste("Step: Getting gene IDs"))

                deGenes_i <- intersect(deGenes(), rownames(gene_ids_all))
                deGenes_i <- gene_ids_all[deGenes_i, 3]
                lende_i <- length(deGenes_i)

                print(kegg_dispdat)
                for (i in 1:nrow(kegg_dispdat)) {
                    incProgress(1/3 + 1/nrow(kegg_dispdat),
                                detail = paste("Step: Calculating "))
                    # cache
                    deGenes_i_ids <- intersect(deGenes_i, kegg_genelists[[i]])
                    leni <- length(deGenes_i_ids)
                    lenkegg <- length(kegg_genelists[[i]])

                    kegg_dispdat[i, 2] <- lenkegg
                    kegg_dispdat[i, 3] <- leni
                    kegg_dispdat[i, 4] <- phyper(leni, lenkegg, nc-1-lenkegg,
                                                 lende_i, lower.tail = F)

                    int_genes <- gene_ids_all[deGenes_i_ids, 1]

                    if (length(int_genes) > 0) {
                        kegg_dispdat[i, 5] <- (paste(int_genes, collapse = ", "))
                    } else {
                        kegg_dispdat[i, 5] <- as.character("0")
                    }
                }

                incProgress(1/3, detail = paste("Step: Formating"))

                kegg_ord <- kegg_dispdat[which(kegg_dispdat[,4] < 0.05),]
                kegg_ord <- kegg_ord[order(kegg_ord[, 4]),]
                kegg_ord[, 4] <- format(kegg_ord[, 4], scientific = T)

                showNotification("KEGG calculation finished")

                #output$downloadKEGG <- downloadHandler(
                #    filename = function() {
                #        paste("KEGG_data", ".csv", sep = "")
                #    },
                #    content = function(file) {
                #        write.csv(kegg_ord, file, row.names = FALSE)
                #    }
                #)
                return(head(kegg_ord, n = 10))
            })
          }, bordered = T)
    }})
}
