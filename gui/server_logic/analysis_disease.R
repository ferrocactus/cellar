load('gui/obj/gene_ids_all')
load('gui/obj/disease_genes')



analysis_disease <- function(input, output, session, deGenes, pipe) {

    dis_names <- names(disease_genes)

    len_dis <- length(dis_names)
    hypergeom_dis <- data.frame(dis_names,
                            integer(length = len_dis),
                            integer(length = len_dis),
                            double(length = len_dis))
    marker_genes_dis <- character(length = nrow(hypergeom_dis))
    hypergeom_dis <- data.frame(hypergeom_dis, marker_genes_dis, stringsAsFactors = F)
    colnames(hypergeom_dis) <- c("Name", "n", "Intersection Length",
                            "pval", "Intersection Genes")

    rownames(gene_ids_all) <- gene_ids_all[,1]

    observe({
    if (length(deGenes()) > 0 && pipe() != 0) {
        ns <- length(pipe()$col_ids)
        s1=as.character(input$subset1)
        s2=as.character(input$subset2)
        tabletitle=paste(s1," (vs. ",s2,")",sep="")
        output$titleD <- renderText(tabletitle)
        output$disease <- DT::renderDataTable({
            withProgress(message = 'Finding Markers', value = 0, {
                incProgress(1/3, detail = paste("Step: Getting gene IDs"))
                deGenes_i_dis <- intersect(deGenes(), rownames(gene_ids_all))

                incProgress(1/3, detail = paste("Step: Calculating"))
                for (i in 1:nrow(hypergeom_dis)) {
                    deGenes_i_ids_dis <- intersect(deGenes_i_dis,
                                               disease_genes[[i]])
                    leni_dis <- length(deGenes_i_ids_dis)
                    lenhs_dis <- length(disease_genes[[i]])

                    hypergeom_dis[i, 1] <- names(disease_genes)[i]
                    hypergeom_dis[i, 2] <- lenhs_dis
                    hypergeom_dis[i, 3] <- leni_dis

                    if (hypergeom_dis[i, 3] == 0){
                        hypergeom_dis[i, 4] <- 1.0
                    } else {
                        hypergeom_dis[i,4]<-phyper(leni_dis, lenhs_dis, ns-1-lenhs_dis,
                                            length(deGenes_i_dis), lower.tail = F)

                        if (length(deGenes_i_ids_dis)>0){
                            hypergeom_dis[i, 5] <- (paste(deGenes_i_ids_dis,
                                                      collapse = ", "))
                        } else {
                            hypergeom_dis[i, 5] <- as.character("0")
                        }
                    }
                }

                ord_dis <- hypergeom_dis[which(hypergeom_dis[,4] < 0.05),]
                ord_dis <- ord_dis[order(ord_dis[, 4]),]
                ord_dis[, 4] <- format(ord_dis[,4], scientific = T)

                showNotification("Disease Markers analysis finished")

                output$downloadDis <- downloadHandler(
                    filename = function() {
                        paste0(pipe()$dataset, "_Disease", "_analysis", ".csv")
                    },
                    content = function(file) {
                        write.csv(ord_dis, file, row.names = FALSE)
                    }
                )
                return(ord_dis)
            })
        }#, bordered = T
        )
    }})
}
