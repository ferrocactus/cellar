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

build_heatmap <- function(adata, heatmap_var) {
    degenes<-py_to_r(get_gene_names_de(adata()))

    if (length(degenes) != 0){
        withProgress(message='Constructing heatmap', {
            label_names=py_to_r(get_labels(adata()))

            cluster_labs=as.character(py_to_r(get_cluster_label_list(adata())))
            incProgress(1 / 3)

            #create heatmap object
            heatmap_dat<-matrix(nrow = length(cluster_labs),ncol = length(degenes))
            #labels for heatmap
            rownames(heatmap_dat)<-cluster_labs
            colnames(heatmap_dat)<-degenes
            scdata_subset<-data.frame(label_names, py_to_r(get_x(adata())))

            colnames(scdata_subset)=c("cluster", py_to_r(get_all_gene_names(adata())))
            #populate the heatmap object
            for (i in 1:length(cluster_labs)) {
                #heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset$clusters==cluster_labs[i]),degenes])
                heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset[,1]==as.double(cluster_labs[i])),degenes])
            }
            incProgress(1 / 3)

            #set the color scale
            scaleRYG <- colorRampPalette(c("blue","white","red"), space = "rgb")(30)
            #plot the heatmap

            heatmap_var(heatmap.2(heatmap_dat, density.info = "none",trace = "none",col = scaleRYG,
                        xlab = "DE Genes",margins = c(9,7),
                        ylab = "Cluster"))
            incProgress(1 / 3)
        })
        return(heatmap_var())
    }
    return(NULL)
}

analysis_body <- function(input, output, session, adata, deGenes, activeDataset) {
    markers_genelists_list <- getMarkerGeneList("src/markers/cell_type_marker.json")
    uploaded_file_flag <- reactiveVal(0)
    heatmap_var = reactiveVal(NULL)

    observe({
        # only run if deGenes have been stored
        req(adata())
        if (length(deGenes()) < 1) return()

        dataset = as.character(activeDataset())
        nc = py_to_r(adata()$shape[1])
        #nc = py_to_r(adata()$n_vars)

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

        output$titleheatmap <- renderText(tabletitle)
        output$heatmap <- renderPlot({
            build_heatmap(adata, heatmap_var)
        }, height=input$heat_height)
    })

    observeEvent(input$markjson, {
        isolate(uploaded_file_flag(uploaded_file_flag() + 1))
    })

    observeEvent(input$heat_height, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$uns, 'de')) return()

        output$heatmap <- renderPlot({
            return(heatmap_var())
        }, height = input$heat_height)
    })

    observeEvent(input$markjson, {
        isolate(uploaded_file_flag(uploaded_file_flag() + 1))
    })
    trigger_threshold <- reactiveVal(FALSE)
    listen_violin <- reactive({
        list(input$color, input$violin_t, input$switcher)
    })

    observeEvent(input$color,{
        trigger_threshold(FALSE)
    })
    observeEvent(input$violin_t,{
        trigger_threshold(TRUE)
    })
<<<<<<< HEAD

    # observeEvent(listen_violin() ,{
    #     if (input$color!='Uncertainty' && input$color != 'Clusters'){
    #         gene_names = py_to_r(get_all_gene_names(adata()))
    #         selected_gene=input$color
    #         i = which(gene_names == (selected_gene))[1]


    #         if (i>0){

    #             gene_data = py_to_r((adata()$X$T[i]))

    #             # index=as.integer(length(gene_data)/10)
    #             # print(index)
    #             # tenp=sort(gene_data)[index]
    #             # print(tenp)
    #             # gene_data=gene_data/tenp

    #             #EPS=0.01
    #             # m = signif(min(gene_data) - EPS, digits=3)
    #             # M = signif(max(gene_data) + EPS, digits=3)
    #             #
    #             # if (isolate(trigger_threshold()) == TRUE) {
    #             #     v1 = isolate(input$violin_t)[1]
    #             #     v2 = isolate(input$violin_t)[2]
    #             # } else {
    #             #     v1 = m
    #             #     v2 = M
    #             # }
    #             v1 = isolate(input$violin_t)[1]
    #             v2 = isolate(input$violin_t)[2]
    #             v1=as.numeric(v1)
    #             v2=as.numeric(v2)
    #             if (v1==4.99 && v2==5.11)
    #                 return()
    #             #index1=which(gene_data %in% gene_data[gene_data>v1])
    #             index1=which(gene_data %in% gene_data[gene_data>v1])# && gene_data %in% gene_data[gene_data>0])
    #             index2=which(gene_data %in% gene_data[gene_data<v2])
    #             index=c()
    #             for (i in 1:length(gene_data)){
    #                 if (i %in% index1 && i %in% index2){
    #                     index=c(index,i)
    #                 }
    #             }
    #             violin_dat0 = data.frame(as.factor(py_to_r(get_labels(adata()))),gene_data)
    #             violin_dat = data.frame(as.factor(py_to_r(get_labels(adata())))[index],gene_data[index])
    #             colnames(violin_dat0)=c("cluster","expression")
    #             colnames(violin_dat)=c("cluster","expression")


    #             setname <- py_to_r(get_label_names(adata()))

    #             setname=setname[index]
    #             lvls=length(levels(as.factor(setname)))
    #             if (length(index)==0 || lvls==length(setname)){
    #                 showNotification("No cell in the thresholds")
    #                 output$violin<-NULL
    #                 output$zeros <-NULL
    #                 return
    #             }

    #             else{

    #                     output$violin <- renderPlotly({
    #                         #ylim1 = boxplot.stats(gene_data)$stats[c(1, 5)]# scale y limits based on ylim1
    #                         ggplot(violin_dat, aes(x=cluster, y=expression, fill=setname))+   geom_violin( )
    #                         #p1 = p0 + coord_cartesian(ylim = ylim1*1.05)

    #                     })
    #                     output$zeros <-renderPlotly({

    #                         data_split = split(violin_dat0,violin_dat0$cluster)
    #                         num_less_zero = matrix(nrow = length(data_split),ncol = 2)
    #                         colnames(num_less_zero) = c("cluster","percentage")
    #                         num_less_zero = data.frame(num_less_zero)
    #                         num_less_zero[,1] = names(data_split)
    #                         for (i in c(1:nrow(num_less_zero))){
    #                             num_less_zero[i,2] = round(sum(data_split[[num_less_zero[i,1]]]$expression<=0)*100/nrow(data_split[[num_less_zero[i,1]]]),3)
    #                         }
    #                         num_less_zero$cluster = as.double(num_less_zero$cluster)
    #                         ggplot(data=num_less_zero, aes(x=cluster, y=percentage)) +
    #                             geom_bar(stat="identity", fill="steelblue")+
    #                             geom_text(aes(label=percentage), vjust=1.6, color="black", size=3.5)+
    #                             #theme_minimal()+
    #                             scale_x_continuous(breaks = seq(0, max(num_less_zero$cluster), 1))+
    #                             ggtitle("0 expression cells")
    #                     })
    #             }
    #             viotitle=paste0("Violin Plot for ",as.character(selected_gene))
    #             output$titleviolin <- renderText(viotitle)
    #         }
    #         else{
    #             showNotification("Gene name not found")
    #         }
    #     }
    #     ## violin
    # })

=======
    
    observeEvent(listen_violin() ,{
        if (input$color!='Uncertainty' && input$color != 'Clusters' && input$switcher=='Violin Plot'){

            v1 = isolate(input$violin_t)[1]
            v2 = isolate(input$violin_t)[2]
            v1=as.numeric(v1)
            v2=as.numeric(v2)
            
            if (v1==4.99 && v2==100)
                return()
            
            progress='Constructing violin plot'
            if (v1==-1 && v2==10)
                progress='Initializing violin plot'            
            withProgress(message=progress, {
            #print(input$switcher)
            
            gene_names = py_to_r(get_all_gene_names(adata()))
            selected_gene=input$color
            i = which(gene_names == (selected_gene))[1]
            gene_data = py_to_r((adata()$X$T[i]))
            lbls=adata()$obs['labels']
            #print(lbls)
            
            incProgress(1 / 3)
            #status=generate_violin(r_to_py(adata()),as.character(input$color),v1,v2)
            status=generate_violin(r_to_py(gene_data),lbls,as.character(input$color),v1,v2)
            #print(status)
            if (status==-1){
                #showNotification("Not enough cells in the thresholds")
                output$violin<-NULL
            }
            else{
                output$violin <- renderImage({
                    list(src = paste0('violin',as.character(input$color),'.png'),
                         contentType = 'image/png',
                         width = 800,
                         height = 600,
                         alt = "This is alternate text")
                }, deleteFile = TRUE)
            }

            
            
            
            
            incProgress(1 / 3)
            if (i>0){

                gene_data = py_to_r((adata()$X$T[i]))
                zero=min(gene_data)

                index1=which(gene_data %in% gene_data[gene_data>v1])# && gene_data %in% gene_data[gene_data>0])
                index2=which(gene_data %in% gene_data[gene_data<v2])
                index=c()
                for (i in 1:length(gene_data)){
                    if (i %in% index1 && i %in% index2){
                        index=c(index,i)
                    }
                }
                violin_dat0 = data.frame(as.factor(py_to_r(get_labels(adata()))),gene_data)
                violin_dat = data.frame(as.factor(py_to_r(get_labels(adata())))[index],gene_data[index])
                colnames(violin_dat0)=c("cluster","expression")
                colnames(violin_dat)=c("cluster","expression")


                setname <- py_to_r(get_label_names(adata()))

                setname=setname[index]
                lvls=length(levels(as.factor(setname)))
                if (length(index)==0 || lvls==length(setname)){
                    #showNotification("No cell in the thresholds")
                    output$violin<-NULL
                    output$zeros <-NULL
                    return
                }

                else{

                        # output$zeros <-renderPlotly({
                        # 
                        #     data_split = split(violin_dat0,violin_dat0$cluster)
                        #     num_less_zero = matrix(nrow = length(data_split),ncol = 2)
                        #     colnames(num_less_zero) = c("cluster","percentage")
                        #     num_less_zero = data.frame(num_less_zero)
                        #     num_less_zero[,1] = names(data_split)
                        #     
                        #     print(zero)
                        #     for (i in c(1:nrow(num_less_zero))){
                        #         num_less_zero[i,2] = round(sum(data_split[[num_less_zero[i,1]]]$expression<=zero)*100/nrow(data_split[[num_less_zero[i,1]]]),5)
                        #     }
                        #     num_less_zero$cluster = as.double(num_less_zero$cluster)
                        #     ggplot(data=num_less_zero, aes(x=cluster, y=percentage)) +
                        #         geom_bar(stat="identity", fill="steelblue")+
                        #         geom_text(aes(label=percentage), vjust=1.6, color="black", size=3.5)+
                        #         #theme_minimal()+
                        #         scale_x_continuous(breaks = seq(0, max(num_less_zero$cluster), 1))+
                        #         ggtitle("0 expression cells")
                        # 
                        # })
                }


                viotitle=paste0("Violin Plot for ",as.character(selected_gene))


                output$titleviolin <- renderText(viotitle)
            }
            })
        }
        else{
            output$violin<-NULL
            output$titleviolin<-NULL
        }
        ## violin
    })
    
>>>>>>> fix_bar
}


