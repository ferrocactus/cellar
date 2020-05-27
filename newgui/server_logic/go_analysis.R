load('gui/obj/Hs.c5')
load('gui/obj/gene_ids_all')

go_analysis <- function(input, output, session, deGenes, pipe) {
    go_categories <- names(Hs.c5)
    len = length(go_categories)
    go_dispdat <- data.frame(go_categories,
                             integer(length = len),
                             integer(length = len),
                             double(length = len),
                             character(length = len),
                             stringsAsFactors = FALSE)
    colnames(go_dispdat) <- c("Name", "n", "Intersection Length",
                              "pval", "Intersection Genes")

    observe({
    # only run if deGenes have been stored
    if (length(deGenes()) > 0) {
        nc = length(pipe()$col_ids)

        output$GeneOntology <- renderTable({
            rownames(gene_ids_all) <- gene_ids_all[, 1]

            withProgress(message = 'Running GO Analysis', detail = NULL,
                         value = 0, {

                incProgress(1/3, detail = paste("Step: Getting gene IDs"))
                deGenes_i <- intersect(deGenes(), rownames(gene_ids_all))
                deGenes_i <- gene_ids_all[deGenes_i, 3]
                lende_i <- length(deGenes_i)

                incProgress(1/3, detail = paste("Step: Calculating "))
                for (i in 1:nrow(go_dispdat)) {
                    # cache
                    deGenes_i_ids <- intersect(deGenes_i, Hs.c5[[i]])
                    lenhs <- length(Hs.c5[[i]])
                    leni = length(deGenes_i_ids)

                    go_dispdat[i, 2] <- lenhs
                    go_dispdat[i, 3] <- leni
                    go_dispdat[i, 4] <- phyper(leni, lenhs, nc - 1 - lenhs,
                                               lende_i, lower.tail = F)
                    rownames(gene_ids_all) <- gene_ids_all[, 3]
                    int_genes <- gene_ids_all[deGenes_i_ids, 1]
                    if (length(int_genes) > 0) {
                        go_dispdat[i, 5] <- (paste(int_genes, collapse=", "))
                    } else {
                        go_dispdat[i, 5] <- as.character("0")
                    }
                }

                incProgress(1/3, detail = paste("Step: Getting Geneontology"))
                go_ord <- go_dispdat[which(go_dispdat[, 4] < 0.05),]
                go_ord <- go_ord[order(go_ord[, 4]),]
                go_ord[, 4] <- format(go_ord[, 4], scientific = T)

                showNotification("GO Analysis finished")

                #output$downloadGO <- downloadHandler(
                #    filename = function() {
                #        paste("GO_data", ".csv", sep = "")
                #    },
                #    content = function(file) {
                #        write.csv(go_ord, file, row.names = FALSE)
                #    }
                #)
                return(head(go_ord, n = 10))
            })
        }, bordered = T)
    }})
}
