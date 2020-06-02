library(rjson)

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

analysis_markers_user <- function(input, output, session, deGenes, pipe) {

    req(input$markjson)
    user_genelists <- getMarkerGeneList(input$markjson$datapath)
    user_gene_names <- names(user_genelists)

    len_user <- length(user_gene_names)
    hypergeom_user <- data.frame(user_gene_names,
                            integer(length = len_user),
                            integer(length = len_user),
                            double(length = len_user))
    marker_genes_user <- character(length = nrow(hypergeom_user))
    hypergeom_user <- data.frame(hypergeom_user, marker_genes_user, stringsAsFactors = F)
    colnames(hypergeom_user) <- c("Name", "n", "Intersection Length",
                            "pval", "Intersection Genes")

    rownames(gene_ids_all) <- gene_ids_all[,1]

    observe({
    if (length(deGenes()) > 0 && pipe() != 0) {
        ns <- length(pipe()$col_ids)

        output$upmarkers <- renderTable({
            withProgress(message = 'Finding Markers', value = 0, {
                incProgress(1/3, detail = paste("Step: Getting gene IDs"))
                deGenes_i_user <- intersect(deGenes(), rownames(gene_ids_all))

                incProgress(1/3, detail = paste("Step: Calculating"))
                for (i in 1:nrow(hypergeom_user)) {
                    deGenes_i_ids_user <- intersect(deGenes_i_user,
                                               user_genelists[[i]])
                    leni_user <- length(deGenes_i_ids_user)
                    lenhs_user <- length(user_genelists[[i]])

                    hypergeom_user[i, 1] <- names(user_genelists)[i]
                    hypergeom_user[i, 2] <- lenhs_user
                    hypergeom_user[i, 3] <- leni_user

                    if (hypergeom_user[i, 3] == 0){
                        hypergeom_user[i, 4] <- 1.0
                    } else {
                        hypergeom_user[i,4]<-phyper(leni_user, lenhs_user, ns-1-lenhs_user,
                                            length(deGenes_i_user), lower.tail = F)

                        if (length(deGenes_i_ids_user)>0){
                            hypergeom_user[i, 5] <- (paste(deGenes_i_ids_user,
                                                      collapse = ", "))
                        } else {
                            hypergeom_user[i, 5] <- as.character("0")
                        }
                    }
                }

                ord_user <- hypergeom_user[which(hypergeom_user[,4] < 0.05),]
                ord_user <- ord_user[order(ord_user[, 4]),]
                ord_user[, 4] <- format(ord_user[,4], scientific = T)

                showNotification("User Markers analysis finished")

                #output$downloadusermks <- downloadHandler(
                #    filename = function() {
                #        paste("User_Markers_data", ".csv", sep = "")
                #    },
                #    content = function(file) {
                #        write.csv(ord_user, file, row.names = FALSE)
                #    }
                #)
                return(ord_user)
            })
        }, bordered = T)
    }})
}
