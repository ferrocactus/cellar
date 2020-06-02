load('gui/obj/gene_ids_all')

getMarkerGeneList <- function(path) {
    marker_genelists <- fromJSON(file = path)
    ##create a new list to store marker genes
    markers_genelists_list <- list()

    for (i in 1:length(marker_genelists)) {
        for (j in 1:length(marker_genelists[[i]])) {
            markers_genelists_list[[paste(
                names(marker_genelists)[i],
                names(marker_genelists[[i]])[j],
                sep = " - "
            )]] <- marker_genelists[[i]][[j]]
        }
    }

    return(markers_genelists_list)
}

analysis_markers <- function(input, output, session, deGenes, pipe) {
    markers_genelists_list <- getMarkerGeneList("markers/cell_type_marker.json")
    marker_list_names <- names(markers_genelists_list)

    len <- length(marker_list_names)
    hypergeom <- data.frame(marker_list_names,
                            integer(length = len),
                            integer(length = len),
                            double(length = len))
    marker_genes <- character(length = nrow(hypergeom))
    hypergeom <- data.frame(hypergeom, marker_genes, stringsAsFactors = F)
    colnames(hypergeom) <- c("Name", "n", "Intersection Length",
                            "pval", "Intersection Genes")

    rownames(gene_ids_all) <- gene_ids_all[,1]

    observe({
    if (length(deGenes()) > 0 && pipe() != 0) {
        ns <- length(pipe()$col_ids)

        output$Markers <- renderTable({
            withProgress(message = 'Finding Markers', value = 0, {
                incProgress(1/3, detail = paste("Step: Getting gene IDs"))
                deGenes_i <- intersect(deGenes(), rownames(gene_ids_all))

                incProgress(1/3, detail = paste("Step: Calculating"))
                for (i in 1:nrow(hypergeom)) {
                    deGenes_i_ids <- intersect(deGenes_i,
                                               markers_genelists_list[[i]])
                    leni <- length(deGenes_i_ids)
                    lenhs <- length(markers_genelists_list[[i]])

                    hypergeom[i, 1] <- names(markers_genelists_list)[i]
                    hypergeom[i, 2] <- lenhs
                    hypergeom[i, 3] <- leni

                    if (hypergeom[i, 3] == 0){
                        hypergeom[i, 4] <- 1.0
                    } else {
                        hypergeom[i,4]<-phyper(leni, lenhs, ns-1-lenhs,
                                            length(deGenes_i), lower.tail = F)

                        if (length(deGenes_i_ids)>0){
                            hypergeom[i, 5] <- (paste(deGenes_i_ids,
                                                      collapse = ", "))
                        } else {
                            hypergeom[i, 5] <- as.character("0")
                        }
                    }
                }

                ord <- hypergeom[which(hypergeom[,4] < 0.05),]
                ord <- ord[order(ord[, 4]),]
                ord[, 4] <- format(ord[,4], scientific = T)

                showNotification("Markers analysis finished")

                #output$downloadMKS <- downloadHandler(
                #    filename = function() {
                #        paste("Markers_data", ".csv", sep = "")
                #    },
                #    content = function(file) {
                #        write.csv(hypergeom_ord, file, row.names = FALSE)
                #    }
                #)
                return(ord)
            })
        }, bordered = T)
    }})
}
