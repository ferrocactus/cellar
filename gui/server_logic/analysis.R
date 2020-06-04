load('gui/obj/Hs.c5')
load('gui/obj/Hs.c2')
load('gui/obj/gene_ids_all')
load('gui/obj/kegg_genelists')
load('gui/obj/keggidtoname')
load('gui/obj/gene_ids_all')

build_table <- function(output, mode, fl, deGenes, nc, alpha, dataset, ns) {
    categories <- names(fl)
    if (mode == 'KEGG')
        categories <- kegg_id_toname[categories]

    len = length(categories)
    dispdat <- data.frame(categories,
                          integer(length = len),
                          integer(length = len),
                          double(length = len),
                          character(length = len),
                          stringsAsFactors = FALSE)
    colnames(dispdat) <- c("Name", "n", "Intersection Length",
                            "pval", "Intersection Genes")

    rownames(gene_ids_all) <- gene_ids_all[, 1]

    withProgress(message = paste('Running', mode, 'Analysis'), value = 0, {

        incProgress(1/3, detail = paste("Step: Getting gene IDs"))
        deGenes_i <- intersect(deGenes(), rownames(gene_ids_all))
        deGenes_i <- gene_ids_all[deGenes_i, 3]
        lende_i <- length(deGenes_i)

        incProgress(1/3, detail = paste("Step: Calculating"))
        for (i in 1:nrow(dispdat)) {
            # cache
            deGenes_i_ids <- intersect(deGenes_i, fl[[i]])
            lenhs <- length(fl[[i]])
            leni = length(deGenes_i_ids)

            dispdat[i, 2] <- lenhs
            dispdat[i, 3] <- leni

            if (dispdat[i, 3] == 0){
                dispdat[i, 4] <- 1.0
            } else {
                dispdat[i,4]<-phyper(leni, lenhs, nc-1-lenhs,
                                      lende_i, lower.tail = F)

                rownames(gene_ids_all) <- gene_ids_all[, 3]
                int_genes <- gene_ids_all[deGenes_i_ids, 1]

                if (length(int_genes)>0){
                    dispdat[i, 5] <- (paste(int_genes,
                                              collapse = ", "))
                } else {
                    dispdat[i, 5] <- as.character("0")
                }
            }
        }

        ord <- dispdat[which(dispdat[, 4] < alpha),]
        ord <- ord[order(ord[, 4]),]
        ord[, 4] <- format(ord[, 4], scientific = T)

        showNotification(paste(mode, "analysis finished"))

        downloadid = paste0("download", mode)

        output[[downloadid]] <- downloadHandler(
            filename = function() {
                paste0(dataset, "_", mode, "_analysis", ".csv")
            },
            content = function(file) {
                write.csv(ord, file, row.names = FALSE)
            }
        )
        return(ord)
    })
}

analysis <- function(input, output, session, deGenes, pipe) {
    observe({
    # only run if deGenes have been stored
    if (length(deGenes()) > 0 && pipe() != 0) {
        nc = length(pipe()$col_ids)

        s1=as.character(isolate(input$subset1))
        s2=as.character(isolate(input$subset2))
        tabletitle=paste(s1," (vs. ",s2,")",sep="")
        output$titleONTO <- renderText(tabletitle)
        output$GOtable <- DT::renderDataTable({
            build_table(output = output, mode = 'GO', fl = Hs.c5,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = pipe()$dataset, ns = session$ns)
        },
        #caption = tabletitle,
        #caption.placement = getOption("xtable.caption.placement", "top"),
        #caption.width = getOption("xtable.caption.width", NULL),
        #bordered = T
        )
        output$titleKEGG <- renderText(tabletitle)
        output$KEGGtable <- DT::renderDataTable({
            build_table(output = output, mode = 'KEGG', fl = kegg_genelists,
                        deGenes = deGenes,
                        nc = nc, alpha = as.numeric(input$mark_alpha),
                        dataset = pipe()$dataset, ns = session$ns)
        }
        #caption = tabletitle,
        #caption.placement = getOption("xtable.caption.placement", "top"),
        #caption.width = getOption("xtable.caption.width", NULL)
        #,bordered = T
        )
        output$titleMSIG <- renderText(tabletitle)
        output$MSIGDBtable <- DT::renderDataTable({
            build_table(output = output, mode = 'MSIGDB', fl = Hs.c2,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = pipe()$dataset, ns = session$ns)
        },
        #caption = tabletitle,
        #caption.placement = getOption("xtable.caption.placement", "top"),
        #caption.width = getOption("xtable.caption.width", NULL),
        #bordered = T
        )
    }})
}
