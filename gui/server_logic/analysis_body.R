load('gui/obj/Hs.c5')
load('gui/obj/Hs.c2')
load('gui/obj/gene_ids_all')
load('gui/obj/kegg_genelists')
load('gui/obj/keggidtoname')
load('gui/obj/disease_genes')

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
        if (mode != 'CellType' && mode != 'UCellType' && mode != 'Disease')
            deGenes_i <- gene_ids_all[deGenes_i, 3]
        lende_i <- length(deGenes_i)

        incProgress(1/3, detail = paste("Step: Calculating"))
        for (i in 1:nrow(dispdat)) {
            # cache
            deGenes_i_ids <- intersect(deGenes_i, fl[[i]])
            lenhs <- length(fl[[i]])
            leni = length(deGenes_i_ids)

            if (mode == 'CellType' || mode == 'UCellType' || mode == 'Disease')
                dispdat[i, 1] <- names(fl)[i]
            dispdat[i, 2] <- lenhs
            dispdat[i, 3] <- leni

            if (dispdat[i, 3] == 0){
                dispdat[i, 4] <- 1.0
            } else {
                dispdat[i, 4] <- phyper(leni, lenhs, nc-1-lenhs,
                                        lende_i, lower.tail = F)

                if (mode != 'CellType' && mode != 'UCellType' && mode != 'Disease') {
                    rownames(gene_ids_all) <- gene_ids_all[, 3]
                    deGenes_i_ids <- gene_ids_all[deGenes_i_ids, 1]
                }

                if (length(deGenes_i_ids)>0){
                    dispdat[i, 5] <- (paste(deGenes_i_ids, collapse = ", "))
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

analysis_body <- function(input, output, session, adata, deGenes) {
    markers_genelists_list <- getMarkerGeneList("markers/cell_type_marker.json")
    uploaded_file_flag <- reactiveVal(0)

    observe({
        # only run if deGenes have been stored
        if (is_active(adata()) == FALSE) return()
        if (length(deGenes()) < 1) return()

        dataset = as.character(adata()$uns[['dataset']])

        nc = py_to_r(adata()$n_vars)

        s1 = as.character(isolate(input$subset1))
        s2 = as.character(isolate(input$subset2))
        tabletitle = paste(s1, " (vs. ", s2, ")")

        output$titleONTO <- renderText(tabletitle)
        output$GOtable <- DT::renderDataTable({
            build_table(output = output, mode = 'GO', fl = Hs.c5,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })

        output$titleKEGG <- renderText(tabletitle)
        output$KEGGtable <- DT::renderDataTable({
            build_table(output = output, mode = 'KEGG', fl = kegg_genelists,
                        deGenes = deGenes,
                        nc = nc, alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })

        output$titleMSIG <- renderText(tabletitle)
        output$MSIGDBtable <- DT::renderDataTable({
            build_table(output = output, mode = 'MSIGDB', fl = Hs.c2,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })

        output$titleCellType <- renderText(tabletitle)
        output$CellTypetable <- DT::renderDataTable({
            build_table(output = output, mode = 'CellType',
                        fl = markers_genelists_list,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })

        output$titleUCellType <- renderText(tabletitle)
        output$UCellTypetable <- DT::renderDataTable({
            if (uploaded_file_flag() < 1) return()
            req(input$markjson)
            user_genelists <- getMarkerGeneList(input$markjson$datapath)
            build_table(output = output, mode = 'UCellType',
                        fl = user_genelists,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })

        output$titleDisease <- renderText(tabletitle)
        output$Diseasetable <- DT::renderDataTable({
            build_table(output = output, mode = 'Disease',
                        fl = disease_genes,
                        deGenes = deGenes, nc = nc,
                        alpha = as.numeric(input$mark_alpha),
                        dataset = dataset, ns = session$ns)
        })
    })

    observeEvent(input$markjson, {
        isolate(uploaded_file_flag(uploaded_file_flag() + 1))
    })
}
