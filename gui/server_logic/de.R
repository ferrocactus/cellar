differential_e <- function(input, output, session, adata, remark, deGenes) {
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
            cellar$de(
                x = adata(),
                subset1 = s1,
                subset2 = s2,
                alpha = input$mark_alpha,
                max_n_genes = input$mark_markers_n,
                correction = input$mark_correction,
                inplace = TRUE)
        })

        remark(remark() + 1)
    })

    observe({
        if (remark() < 1) return()
        isolate(remark(0))

        # Switch to DE table
        #reset_analysis_tabs(output)
        updateTabsetPanel(session, "switcher", selected = "DE")
        deGenes(c())

        ns <- session$ns

        #reset_analysis_tabs(output)
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
                py_to_r(get_gene_logFC_de(adata())),
                py_to_r(get_gene_pvals_de(adata())),
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
            table <- table[order(table[, 2], decreasing = TRUE),]
            colnames(table) <- c(
                "DE genes", "logFC", "adj.pval", "Show Expression Level")
            table[, 3] <- format(table[, 3], scientific = T)

            deGenes(table[, 1])

            output$downloadDE <- downloadHandler(
                filename = function() {
                    paste0("Cellar", "_DE_genes", ".csv")
                },
                content = function(file) {
                    write.csv(table, file, row.names = FALSE)
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

        # observeEvent(input$heat_height,{
        #     degenes<-pipe()$markers[['0']][['outp_names']]
        #     if (length(degenes) != 0){
        #         output$titleheatmap <- renderText(tabletitle)
        #         #label_names = pipe()$get_label_names()
        #         label_names=pipe()$labels
        #         #get the cluster names
        #         cluster_labs<-as.character(levels(as.factor(pipe()$labels)))
        #         #get the de genes
        #         #create heatmap object
        #         heatmap_dat<-matrix(nrow = length(cluster_labs),ncol = length(degenes))
        #         #labels for heatmap
        #         rownames(heatmap_dat)<-cluster_labs
        #         colnames(heatmap_dat)<-degenes
        #         scdata_subset<-data.frame(pipe()$labels, pipe()$x)
        #         colnames(scdata_subset)=c("cluster", pipe()$col_names)
        #         #populate the heatmap object
        #         for (i in 1:length(cluster_labs)) {
        #             #heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset$clusters==cluster_labs[i]),degenes])
        #             heatmap_dat[cluster_labs[i],]<-colMeans(scdata_subset[which(scdata_subset[,1]==as.double(cluster_labs[i])),degenes])
        #         }
        #         output$heatmap <- renderPlot({
        #             #set the color scale
        #             scaleRYG <- colorRampPalette(c("blue","white","red"), space = "rgb")(30)
        #             #plot the heatmap

        #             heatmap.2(heatmap_dat,density.info = "none",trace = "none",col = scaleRYG,
        #                       xlab = "DE Genes",margins = c(9,7),
        #                       ylab = "Cluster")
        #         },height=input$heat_height
        #         )
        #     }
        # })

        ## end of heatmap
    })
}

reset_analysis_tabs <- function(output) {
    # Reset analysis tabs and switch to DE
    # By default the tabs will be erased only when we switch to them
    # so we first disable suspendWhenHidden, set them to null, and enable
    # it again
    outputOptions(output, "DEtable", suspendWhenHidden = FALSE)
    outputOptions(output, "KEGGtable", suspendWhenHidden = FALSE)
    outputOptions(output, "GOtable", suspendWhenHidden = FALSE)
    outputOptions(output, "CellTypetable", suspendWhenHidden = FALSE)
    outputOptions(output, "MSIGDBtable", suspendWhenHidden = FALSE)
    outputOptions(output, "UCellTypetable", suspendWhenHidden = FALSE)
    outputOptions(output, "Diseasetable", suspendWhenHidden = FALSE)
    outputOptions(output, "heatmap", suspendWhenHidden = FALSE)
    output$DEtable = NULL
    output$GOtable = NULL
    output$KEGGtable = NULL
    output$CellTypetable = NULL
    output$UCellTypetable = NULL
    output$MSIGDBtable = NULL
    output$Diseasetable = NULL
    output$heatmap = NULL
    outputOptions(output, "DEtable", suspendWhenHidden = TRUE)
    outputOptions(output, "KEGGtable", suspendWhenHidden = TRUE)
    outputOptions(output, "GOtable", suspendWhenHidden = TRUE)
    outputOptions(output, "CellTypetable", suspendWhenHidden = TRUE)
    outputOptions(output, "MSIGDBtable", suspendWhenHidden = TRUE)
    outputOptions(output, "UCellTypetable", suspendWhenHidden = TRUE)
    outputOptions(output, "Diseasetable", suspendWhenHidden = TRUE)
    outputOptions(output, "heatmap", suspendWhenHidden = TRUE)
}