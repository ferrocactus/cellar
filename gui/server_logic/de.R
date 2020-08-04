differential_e <- function(input, output, session, adata, remark, deGenes,
                           activeDataset) {
    observeEvent(input$getdegenes, {
        if (is_active(adata()) == FALSE) return()
        if (!py_has_attr(adata()$obs, 'labels')) return()

        s1 = as.character(input$subset1)
        s2 = as.character(input$subset2)


        if (s1 == 'None' && s2 == 'None') {
            showNotification("Please select a subset to analyze.")
            return()
        }
        if (s1 == s2) {
            showNotification("The selected subsets should be different.")
            return()
        }

        if (s1 == 'None') {
            s1 = s2
            s2 = NULL
        } else if (s2 == 'None') {
            s2 = NULL
        }

        withProgress(message = "Please Wait", value = 0, {
            n <- 2

            incProgress(1 / n, detail = "Finding DE Genes")
            if (input$color=='Uncertainty'){
                msg <- cellar$safe(cellar$de,
                                   x = adata(),
                                   subset1 = s1,
                                   subset2 = s2,
                                   method = input$test_method,
                                   alpha = input$mark_alpha,
                                   max_n_genes = input$mark_markers_n,
                                   correction = input$mark_correction,
                                   uncertain=TRUE,
                                   inplace = TRUE)
            }
            else{
                msg <- cellar$safe(cellar$de,
                                   x = adata(),
                                   subset1 = s1,
                                   subset2 = s2,
                                   method = input$test_method,
                                   alpha = input$mark_alpha,
                                   max_n_genes = input$mark_markers_n,
                                   correction = input$mark_correction,
                                   inplace = TRUE)
            }


            if (msg != 'good') {
                showNotification(py_to_r(msg))
                return()
            }
        })

        remark(remark() + 1)
    })

    observe({
        if (remark() < 1) return()
        isolate(remark(0))
        if (has_key(adata(), 'uns', 'de') == FALSE) return()

        # Switch to DE table
        reset_analysis_tabs(output)
        updateTabsetPanel(session, "switcher", selected = "DE")
        isolate(deGenes(c()))

        ns <- session$ns

        reset_analysis_tabs(output)
        s1 = as.character(input$subset1)
        s2 = as.character(input$subset2)
        tabletitle = paste0(s1, " (vs. ", s2, ")")
        output$titleDE <- renderText(tabletitle)

        shinyInput <- function(FUN, len, id, ns, ...) {  ## for inputting buttons
            inputs <- character(len)
            for (i in seq_len(len)) {
                inputs[i] <- as.character(FUN(paste0(id, i), ...))
            }
            inputs
        }

        de_gene_names = py_to_r(get_gene_names_de(adata()))

        output$DEtable <- DT::renderDataTable({
            table <- data.frame(
                de_gene_names,
                py_to_r(get_de_table(adata())),
                Actions = shinyInput(
                    actionButton,
                    length(de_gene_names),
                    'button_',
                    label = "Show Expression",
                    onclick = sprintf("Shiny.onInputChange('%s', this.id)",
                                      ns("select_button")),
                ),
                check.names = TRUE, # remove duplicates etc
                stringsAsFactors = FALSE
            )
            table <- table[order(table$'log2fc', decreasing = TRUE),]
            #colnames(table) <- c(
            #    "DE genes", "logFC", "adj.pval", "Show Expression Level")
            table <- format(table, scientific=T)

            isolate(deGenes(table[, 1]))

            dataset = as.character(activeDataset())
            output$downloadDE <- downloadHandler(
                filename = function() {
                    paste0(dataset, "_DE_genes", ".csv")
                },
                content = function(file) {
                    write.csv(table[,1:ncol(table)-1], file, row.names = FALSE)
                }
            )

            return(table)
        }, server = FALSE, escape = FALSE, selection = 'none')

        o <- observeEvent(input$select_button, {
            if (identical(input$select_button, 'reset')) {
                return()
            }
            selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
            if (length(de_gene_names) >= selectedRow) {
                selected_gene = de_gene_names[[selectedRow]]
                updateSelectInput(
                    session,
                    'color',
                    selected = selected_gene
                )
            }
        })

        #################### heat map
        ## heatmap

        observeEvent(input$heat_height,{
            if (is_active(adata()) == FALSE) return()
            if (!py_has_attr(adata()$obs, 'labels')) return()
            degenes<-py_to_r(get_gene_names_de(adata()))

            if (length(degenes) != 0){
                output$titleheatmap <- renderText(tabletitle)
                #label_names = pipe()$get_label_names()
                label_names=py_to_r(get_labels(adata()))

                cluster_labs=as.character(py_to_r(get_cluster_label_list(adata())))
                #get the cluster names



                #get the de genes
                #create heatmap object
                heatmap_dat<-matrix(nrow = length(cluster_labs),ncol = length(degenes))
                #labels for heatmap
                rownames(heatmap_dat)<-cluster_labs
                colnames(heatmap_dat)<-degenes
                scdata_subset<-data.frame(label_names, py_to_r(adata()$X))

                colnames(scdata_subset)=c("cluster", py_to_r(get_all_gene_names(adata())))
                #populate the heatmap object
                for (i in 1:length(cluster_labs)) {
                    #heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset$clusters==cluster_labs[i]),degenes])
                    heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset[,1]==as.double(cluster_labs[i])),degenes])
                }
                output$heatmap <- renderPlot({
                    #set the color scale
                    scaleRYG <- colorRampPalette(c("blue","white","red"), space = "rgb")(30)
                    #plot the heatmap

                    heatmap.2(heatmap_dat,density.info = "none",trace = "none",col = scaleRYG,
                              xlab = "DE Genes",margins = c(9,7),
                              ylab = "Cluster")
                },height=input$heat_height
                )
            }
        })

        ## end of heatmap
    })
}

reset_analysis_tabs <- function(output) {
    output$DEtable = NULL
    output$GOtable = NULL
    output$KEGGtable = NULL
    output$CellTypetable = NULL
    output$UCellTypetable = NULL
    output$MSIGDBtable = NULL
    output$Diseasetable = NULL
    output$heatmap = NULL
}
