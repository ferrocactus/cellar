library(reticulate)
library(shiny)
library(plotly)
library(ggplot2)
library(stringr)
library(limma)
library(GO.db)
library(rjson)

source_python('__init__.py')
load('gui/obj/Hs.c2')
load('gui/obj/Go_kegg_lists')
load('gui/obj/cell_ontology')
load('gui/obj/keggidtoname')
load('gui/obj/gene_ids_all')

# Load utility functions
source("gui/functions.R") # getPage, intersect, writeDataset, getHypergeom
                          # dataModal
source("gui/pipe_functions.R") # runPipe

################################################################# server
server <- shinyServer(function(input, output, session) {
  #SESSION-WISE VARIABLES
  env=environment()

  new_labels=NULL
  assign("new_labels",NULL,envir = env)


  # count of selected sets
  setcount=0
  assign("setcount", 0, envir = env)


  degenenames=NULL
  assign("degenenames",NULL,envir = env)
  firstflag=1
  assign("firstflag",1,envir = env)
  clusters=0
  assign("clusters",0,envir = env)
  debuttons=NULL
  assign("debuttons",NULL,envir=env)
  intersectbuttons=NULL
  switcher=NULL
  assign("intersectbuttons",NULL,envir = env)
  assign("switcher",NULL,envir = env)
  sst=NULL
  assign("sst",NULL, envir = env)
  ss=FALSE
  assign("ss",FALSE, envir=env)
  dataset = "default"
  assign("dataset", "default", envir=env)
  #END OF SESSION-WISE VARIABLES

  ###########################################################
  # UI Logic
  observeEvent(input$dim_method, {
      if(input$dim_method == "PCA"){
          shinyjs::enable(selector = "[type=radio][value=pca_auto]")
      } else {
          shinyjs::disable(selector = "[type=radio][value=pca_auto]")
          updateRadioButtons(
              session,
              "dim_options",
              selected = "pca_manual"
          )
      }
  })

  observeEvent(input$dim_options, {
      if(input$dim_options == "pca_manual"){
          shinyjs::enable('dim_n_components')
      } else {
          shinyjs::disable('dim_n_components')
      }
  })

  observeEvent(input$clu_method, {
      if(input$clu_method %in% options$clu_no_n_clusters){
          shinyjs::disable('clu_n_clusters')
          shinyjs::disable('eval_method')
      } else {
          shinyjs::enable('clu_n_clusters')
          shinyjs::enable('eval_method')
      }
  })

  observeEvent(input$ensemble_checkbox, {
      if('all' %in% str_split(input$ensemble_checkbox, " ")){
          for (x in options$clu_ensemble) {
              if (x != 'all') {
                  shinyjs::disable(
                      selector = paste("#ensemble_checkbox input[value='",
                                  x, "']", sep="")
                  )
              }
          }
      } else {
          for (x in options$clu_ensemble) {
              if (x != 'all') {
                  shinyjs::enable(
                      selector=paste("#ensemble_checkbox input[value='",x, "']",
                                       sep="")
                  )
              }
          }
      }
  })
  ###########################################################

  # Upload dataset
  # when a dataset is uploaded, write it into the "datasets" folder
  # users can choose it when run another configuration

  observeEvent(input$file1, {
    req(input$file1)

    tryCatch({
        # Dialog window
        showModal(dataModal(code = 0))

    }, error = function(e) {
        stop(safeError(e))
    })
  })

  o <- observeEvent(input$okdataset, {
      # Check if dataset exists
      fname <- input$dataset_fname
      print(fname)

      if (fname == "") {
          showModal(dataModal(code = 2))
      } else if (datasetExists(fname, path="datasets/user_uploaded")) {
          showModal(dataModal(code = 1))
      } else {
          removeModal()
          writeDataset(input$file1$datapath, fname)
          updateSelectInput(
              session = session,
              inputId = "uploaded_dataset",
              label = "Choose a dataset:",
              choices = list.files("datasets/user_uploaded"),
              showNotification("Dataset uploaded")
          )
      }
      #removeUI(selector = "#dataset_fname")
      o$destroy
  })


  # rerun the app if "run with new configuration button" pressed
  # (this can avoid observing previous events)
  observeEvent(input$runconfigbtn, {

    assign("ss", FALSE, envir = env)

    if (input$folder=="user_uploaded"){
      dataset=as.character(input$uploaded_dataset)
      print(dataset)
      showNotification(paste("Dataset: ",dataset,sep=""))
      pipe <- Pipeline(x = dataset, dataset_source=input$folder)
    }
    else{
      dataset=as.character(input$hubmap_dataset)
      showNotification(paste("Dataset: ",dataset,sep=""))
      pipe <- Pipeline(x = dataset, dataset_source=input$folder)
    }

  #  print(dataset)


    tryCatch({
      df <- isolate(runPipe(pipe, input))
    }, error = function(e) {
      df <- "An error occurred."
    }
    )

    if (is.character(df) & length(df) == 1) {
      showNotification(df)
    } else {
      #assign("df", isolate(runPipe(pipe, input)), envir = env)
      ################################### RUN WITH CURRENT CONFIG

      ################################# DISABLE PREVIOUS EVENTS
      #disable observing buttons
      if (length(sst)>0){
        for (i in 1:length(sst)){
          sst[[i]]$destroy()
        }
      }
      if (length(debuttons)>0){
        for (i in 1:length(debuttons)){
          debuttons[[i]]$destroy()
          #removeUI(selector=paste0("div:has(> #",degenenames[[i]],")"))
          #showNotification(degenenames[[i]])
        }
        #assign("degenenames",NULL,envir = env)
        assign("debuttons",NULL,envir=env)
      }
      if (length(intersectbuttons)>0){
        for (i in 1:length(intersectbuttons)){
          intersectbuttons[[i]]$destroy()
        }
        assign("intersectbuttons",NULL,envir = env)
      }
      if (length(switcher)>0){
        for (i in 1:length(switcher)){
          switcher[[i]]$destroy()
        }
        assign("switcher",NULL,envir = env)
      }

      #delete previous UIs
      if (firstflag==0){
        for (i in 1:clusters){
          removeTab(inputId="switcher", target=as.character(i-1))
        }
        output$DEbuttons<-NULL
        output$plot<-NULL
        output$brush<-NULL
        output$genes<-NULL
        #removeUI(selector="#tabset")
        # removeUI(selector="#DEsubsets")
        # insertUI(
        #   selector = "#genecard",
        #   where = "beforeBegin",
        #   ui = actionButton("getdegenes", "Get DE genes", class="sidebtn")
        # )
        #
        # insertUI(
        #   selector = paste0("#","placeholder"),
        #   where = "beforeBegin",
        #   ui = tabsetPanel(
        #     type = "tabs",
        #     id = "tabset",
        #     tabPanel(
        #       "Main Plot",
        #       h3(textOutput("caption")),
        #       plotlyOutput("plot"),
        #       uiOutput("deinfo"),
        #       verbatimTextOutput("genes"),
        #       uiOutput("DEbuttons")
        #     ),
        #     tabPanel(
        #       "Updated Plot",
        #       verbatimTextOutput("brush"),
        #       plotlyOutput("Plot2")
        #     ),
        #     tabPanel(
        #       "Gene Ontology",
        #       verbatimTextOutput("GeneOntology")
        #     ),
        #     tabPanel(
        #       "KEGG",
        #       verbatimTextOutput("KEGG")
        #     ),
        #     tabPanel(
        #       "Markers Intersect",
        #       verbatimTextOutput("Markers")
        #     ),
        #     tabPanel(
        #       "MSigDB C2",
        #       verbatimTextOutput("Msigdb")
        #     )
        #   )
        # )


        #removeUI(selector="div:has(> #DEbuttons)")

      }

      assign("firstflag",0,envir = env)
      ############################ END OF DISABLing PREVIOUS EVENTS
      observeEvent(input$ssclurun, {
        if(exists("updated_new_labels")) {
          tryCatch({
            df <- runSSClu(pipe, updated_new_labels, input)
          }, error = function(e) {
            df <- "An error occurred."
          })

          if (is.character(df) & length(df) == 1) {
            showNotification(df)
          } else {
            markers <- pipe$markers
            #get the gene expression data
            expr_data = data.frame(
              matrix(pipe$x, ncol = length(pipe$col_ids),
                     dimnames = list(1:length(pipe$x[,1]), pipe$col_ids))
            )
            expr_data$cluster = df$y
            ##read in the marker from JSON
            hypergeom <- getHypergeom("markers/cell_type_marker.json")
            markers_genelists_list <- getMarkerGeneList(
              "markers/cell_type_marker.json")


            newlabs <- df[, 3]
            names(newlabs) <- rownames(df)
            labeldats <- levels(as.factor(df[,3]))

            output$plot <- renderPlotly({
              #factorize cluster labels (discrete instead of continuous)
              if (input$color == "cluster") {
                plotcols = as.factor(expr_data[[input$color]])
              } else {
                plotcols = expr_data[[input$color]]
              }
              plot_ly(
                df, x = df$x1, y = df$x2,
                text = ~paste("label: ", as.factor(df$y)),
                color = plotcols,
                key = row.names(df),
                type = 'scatter',
                mode = 'markers'
              ) %>% layout(dragmode = "lasso",
                           title = paste("Value of ", input$color, sep=""))
            })

            ### updated plot
            observeEvent(input$labelupd, {
              #d <- event_data("plotly_selected")
              keys<-s1
              keysNotification(as.character(input$newlabels))

              newlabs[keys]<<-as.character(input$newlabels)


              assign("updated_new_labels", newlabs, envir = env)
              output$Plot2 <- renderPlotly({
                plot_ly(
                  df, x = df[, 1], y = df[, 2],
                  text = ~paste("label: ", as.factor(newlabs)),
                  color = as.factor(newlabs),
                  type = 'scatter',
                  mode = 'markers'
                ) %>% layout(dragmode = "lasso",
                             title = paste("Value of ", input$labelupd, sep=""))
              })
            })
          }
        }
      })
      ####################### End of clearing previous events



      # the dictionary includes information of each cluster
      # (including DE genes and intersections)
      markers <- pipe$markers
      #get the gene expression data
      expr_data = data.frame(
        matrix(pipe$x, ncol = length(pipe$col_ids),
               dimnames = list(1:length(pipe$x[,1]), pipe$col_ids))
      )
      expr_data$cluster = df$y
      ##read in the marker from JSON
      hypergeom <- getHypergeom("markers/cell_type_marker.json")
      markers_genelists_list <- getMarkerGeneList("markers/cell_type_marker.json")
      #Requires HS.c2 to be loaded.
      msigdb_categories <- names(Hs.c2)
      msigdb_pvals <- double(length = length(msigdb_categories))
      msig_dispdat <- data.frame(msigdb_categories, msigdb_pvals)
      newlabs <- df[, 3]
      assign("newlabs",df[, 3],envir = env)
      names(newlabs) <- rownames(df)

      msigdb_categories <- names(Hs.c2)
      msigdb_n <-integer(length = length(msigdb_categories))
      msigdb_nde<-integer(length = length(msigdb_categories))
      msigdb_pvals <- double(length = length(msigdb_categories))
      msigdb_genes<-character(length = length(msigdb_categories))
      msig_dispdat <- data.frame(msigdb_categories, msigdb_n)
      msig_dispdat<-data.frame(msig_dispdat,msigdb_nde)
      msig_dispdat<-data.frame(msig_dispdat,msigdb_pvals)
      msig_dispdat<-data.frame(msig_dispdat,msigdb_genes,stringsAsFactors=F)
      colnames(msig_dispdat)<-c("NAME","N","Intersect Length","P Value","Intersect genes")

      go_categories <- names(Hs.c5)
      go_n<-integer(length = length(go_categories))
      go_nde<-integer(length = length(go_categories))
      go_pvals <- double(length = length(go_categories))
      go_genes<- character(length = length(go_categories))
      go_dispdat <- data.frame(go_categories, go_n)
      go_dispdat<-data.frame(go_dispdat,go_nde)
      go_dispdat<-data.frame(go_dispdat,go_pvals)
      go_dispdat<-data.frame(go_dispdat,go_genes,stringsAsFactors = F)
      colnames(go_dispdat)<-c("NAME","N","Intersect Length","P Value","Intersect Genes")

      kegg_categories <- kegg_id_toname[names(kegg_genelists)]
      kegg_n<-integer(length = length(kegg_categories))
      kegg_nde<-integer(length = length(kegg_categories))
      kegg_pvals <- double(length = length(kegg_categories))
      kegg_genes<-character(length = length(kegg_categories))
      kegg_dispdat <- data.frame(kegg_categories, kegg_n)
      kegg_dispdat<-data.frame(kegg_dispdat,kegg_nde)
      kegg_dispdat<-data.frame(kegg_dispdat,kegg_pvals)
      kegg_dispdat<-data.frame(kegg_dispdat,kegg_genes,stringsAsFactors = F)
      colnames(kegg_dispdat)<-c("NAME","N","Intersect Length","P Value","Intersect Genes")

      hypergeom <- getHypergeom("markers/cell_type_marker.json")
      marker_genes<-character(length = nrow(hypergeom))
      hypergeom<-data.frame(hypergeom,marker_genes,stringsAsFactors = F)
      colnames(hypergeom)<-c("NAME","N","Intersect Length","P Value","Intersect Genes")

      markers_genelists_list <- getMarkerGeneList(
        "markers/cell_type_marker.json")

      cell_ontology_names<-paste(cell_ont_full$id,cell_ont_full$name,sep = " ")
      labeldats<-c(levels(as.factor(expr_data[,length(expr_data)])),cell_ontology_names)

      ############################################################ End of pipeline data processing


      ############################## start observing UI events





      # initialize subsets, store clusters
      assign("set_name",c("None"),envir = env)
      assign("sets",c(NA),envir = env)
      setcount=0
      for (i in 1:length(markers)){
        #names(subsets)[setcount+1]<-paste("cluster_",as.character(i-1),sep="")
        assign("set_name",c(set_name,(paste("cluster_",as.character(i-1),sep=""))),envir=env)
        keys<-which(expr_data$cluster==(i-1))
        assign("sets",c(sets,list(keys)),envir=env)
        #subsets[[as.character(paste("cluster_",as.character(i-1),sep=""))]]<-keys
        assign("setcount", setcount+1, envir = env)
      }


      #list2env(subsets,env)

      # store selected subsets
      observeEvent(input$store_lasso,{
        d <- event_data("plotly_selected")

        keys=as.numeric(d$key)
        cell_count=length(d$key)
        #showNotification(as.character(input$newsubset),duration=NULL)
        if (identical(d$key,NULL)==TRUE){
          return()
        }
        if (identical(input$newsubset,NULL)==FALSE){
          #names(subsets)[setcount+1]<-as.character(input$newsubset)
          #subsets[[as.character(input$newsubset)]]<-keys
          #subsets[setcount+1]<-keys
          sets<<-c(sets,list(keys))
          ###assign("set_name",c(set_name,as.character(input$newsubset)),envir=env)

          set_name<<-c(set_name,as.character(input$newsubset))

          #showNotification(as.character(subsets[setcount+1]),duration=NULL)
          #showNotification(names(subsets)[setcount+1],duration=NULL)
        }
        else{
          #names(subsets)[setcount+1]<-paste("set",as.character(setcount+1),sep="")
          sets<<-c(sets,list(keys))
          set_name<<-c(set_name,paste("set",as.character(setcount+1),sep=""))
          #set_name=c(set_name,paste("set",as.character(setcount+1),sep=""))
        }
        #list2env(subsets,env)
        showNotification(paste(as.character(cell_count)," cells stored",sep=""))
        assign("setcount", setcount+1, envir = env)

        # select subset1
        updateSelectInput(session = session,
                          inputId = "subset1",
                          label = "Choose Subset 1",
                          choices = set_name,
                          selected = NULL)
        # select subset2
        updateSelectInput(session = session,
                          inputId = "subset2",
                          label = "Choose Subset 2",
                          choices = set_name,
                          selected = NULL)
        #showNotification(as.character(setcount),duration=NULL)
      })


      # select subset1
      updateSelectInput(session = session,
                        inputId = "subset1",
                        label = "Choose Subset 1",
                        choices = set_name,
                        selected = NULL)
      # select subset2
      updateSelectInput(session = session,
                        inputId = "subset2",
                        label = "Choose Subset 2",
                        choices = set_name,
                        selected = NULL)



      # select color value
      updateSelectInput(session = session,
                        inputId = "color",
                        label = "Select colour value:",
                        choices = c("cluster", names(expr_data)),
                        selected = NULL)

      # change label

      updateSelectInput(session = session,
                        inputId = "newlabels",
                        label = "Select label",
                        choices = labeldats,
                        selected = NULL)

      # input new label
      observe({
        x <- input$newlabelbox
        #Can use character(0) to remove all choices
        #labeldats<-reactive({value = labchoices})
        observeEvent(input$labeladd, {
          labeldats <<- union(labeldats, x)
          # Can also set the label and select items
          updateSelectInput(session,
                            "newlabels",
                            label = paste("Select input label", length(x)),
                            choices = labeldats)
        })
      })


      ##### gene card
      observeEvent(input$search, {
        if (input$searchgene %in% names(expr_data)) {
          getPage(input$searchgene) # function in functions.R
        }
        else {
          showNotification("Gene name does not exist.")
        }
      })

      ########################################## PLOTTING
      output$plot <- renderPlotly({
        #factorize cluster labels (discrete instead of continuous)
        if (input$color == "cluster") {
          plotcols = as.factor(expr_data[[input$color]])
        } else {
          plotcols = expr_data[[input$color]]
        }
        plot_ly(
          df, x = df$x1, y = df$x2,
          #width=1040,height=780,
          text = ~paste("label: ", as.factor(df$y)),
          color = plotcols,
          key = row.names(df),
          type = 'scatter',
          mode = 'markers'
        ) %>% layout(dragmode = "lasso",
                     title = paste("Value of ", input$color, sep=""))
      })


      ################################# storing selected cells
      assign("s1", NULL, envir = env)  ##
      assign("s2", NULL, envir = env)  ##
      cell_count=0
      cell_count2=0
       assign("cell_count", 0, envir = env)
       assign("cell_cuont2", 0, envir = env)




      ### updated plot


      observeEvent(input$labelupd, {
        #d <- event_data("plotly_selected")
        keys<-subsets[[as.character(input$subset1)]]
        showNotification(as.character(input$newlabels))
        newlabs[keys]<<-as.character(input$newlabels)
        showNotification("Updating labels")

        assign("updated_new_labels", newlabs, envir = env)
        output$Plot2 <- renderPlotly({
          plot_ly(
            df, x = df[, 1], y = df[, 2],
            text = ~paste("label: ", as.factor(newlabs)),
            color = as.factor(newlabs),
            type = 'scatter',
            mode = 'markers'
          ) %>% layout(dragmode = "lasso",
                       title = paste("Value of ", input$labelupd, sep=""))
        })
      })


      ########################################################################### DE GENE IMPLEMENTATION
      ### update select input for cluster selection
      # updateSelectInput(session = session,
      #                   inputId = "cluforanalysis",
      #                   label = "Select clusters",
      #                   choices = levels(as.factor(expr_data$cluster)),
      #                   selected = NULL)
      #
      #
      #
      #
       scdata_subset=expr_data

      #if (length(intersect(rownames(scdata_subset),gene_ids_all[,1]))>0){
       rownames(gene_ids_all)<-gene_ids_all[,1]
      #} else if (length(intersect(rownames(scdata_subset),gene_ids_all[,2]))>0){
      #  rownames(gene_ids_all)<-gene_ids_all[,2]
      #} else if (length(intersect(rownames(scdata_subset),gene_ids_all[,3]))>0){
      #  rownames(gene_ids_all)<-gene_ids_all[,3]
      #} else {
      #  showNotification("Invalid gene ids, functional analysis cannot be performed")
      #}
      observeEvent(input$getdegenes,{
        selecteddat=NULL
        idx1=which(set_name==as.character(input$subset1))
        idx2=which(set_name==as.character(input$subset2))

        assign("s1", sets[[idx1]], envir = env)
        assign("s2", sets[[idx2]], envir = env)

        cell_count=length(s1)
        cell_count2=length(s2)

        if (as.character(input$subset2)=="None"){   ## one against the rest (can be selected by lasso or cluster selection)
          keys=s1
          mar=pipe$get_markers_subset(
                keys-1, alpha=input$mark_alpha, correction=input$mark_correction,
                markers_n=input$nogenes)
          #['indices', 'pvals', 'diffs', 'inp_names', 'outp_names', 'lvl1_type',
          #'lvl1_sv', 'lvl1_intersec', 'lvl1_total', 'lvl1_all', 'lvl2_type',
          #' 'lvl2_sv', 'lvl2_intersec', 'lvl2_total', 'lvl2_all']
          pval=mar[['0']][['pvals']]
          logFC=mar[['0']][['diffs']]
          gene_names=mar[['0']][['outp_names']]


          #selecteddat<-scdata_subset[as.numeric(keys),1:ncol(scdata_subset)-1]
          #restdat<-scdata_subset[-as.numeric(keys),1:ncol(scdata_subset)-1]
          showNotification("DE genes of subset1 against the rest")

        }
        else if(as.character(input$subset1)=="None"){
          keys=s2
          mar=pipe$get_markers_subset(
                keys-1, alpha=input$mark_alpha, correction=input$mark_correction,
                markers_n=input$nogenes)
          pval=mar[['0']][['pvals']]
          logFC=mar[['0']][['diffs']]
          gene_names=mar[['0']][['outp_names']]
          #selecteddat<-scdata_subset[as.numeric(keys),1:ncol(scdata_subset)-1]
          #restdat<-scdata_subset[-as.numeric(keys),1:ncol(scdata_subset)-1]
          showNotification("DE genes of subset2 against the rest")

        }
        else{                 ## 2 subsets
          keys=s1
          keys2=s2
          mar=pipe$get_markers_subset(
                keys-1, keys2-1, alpha=input$mark_alpha,
                correction=input$mark_correction,
                markers_n=input$nogenes)


          ### calculating DE
          pval=mar[['0']][['pvals']]
          logFC=mar[['0']][['diffs']]
          gene_names=mar[['0']][['outp_names']]

          pval2=mar[['1']][['pvals']]
          logFC2=mar[['1']][['diffs']]
          gene_names2=mar[['1']][['outp_names']]



          #selecteddat<-scdata_subset[as.numeric(keys),1:ncol(scdata_subset)-1]
          #restdat<-scdata_subset[as.numeric(keys2),1:ncol(scdata_subset)-1]
          showNotification("DE genes of 2 subsets ")
        }

        degenes_table<-data.frame(gene_names,logFC,stringsAsFactors=F)
        degenes_table<-data.frame(degenes_table,pval,stringsAsFactors=F)
        degenes_table_ord<-degenes_table[order(degenes_table[,2],decreasing=T),]
        colnames(degenes_table_ord)<-c("DE genes","logFC","adj p.value")

        ############# destroy previous DE gene buttons
        if (length(debuttons)>0){
          for (i in 1:length(debuttons)){
            debuttons[[i]]$destroy()
            #debuttons[[1]]<-NULL
          }
        }
        assign("debuttons",NULL,envir=env)  # disable previous buttons
        ###end of dealing with previous buttons


        #### start calculating DE genes
        #### and realizing related functions
        output$genes <- renderTable({



          #####old DE gene implementations
          # withProgress(message = 'calculating DE genes',detail=NULL, value = 0, {
          #   #exp_genes_mean<-colSums(exp_genes)/nrow(exp_genes)
          #
          #   incProgress(1/5, detail = paste("Step: Fetching 2 sets"))
          #   labelsdat<-as.factor(c(rep("selected",nrow(selecteddat)),rep("notselected",nrow(restdat))))
          #   alldat<-rbind(selecteddat,restdat)
          #   alldat<-data.frame(alldat,labelsdat)
          #   modmat<-model.matrix(~labelsdat,data = alldat)
          #   incProgress(2/5, detail = paste("Step: Processing data"))
          #   t=t(alldat[,1:ncol(alldat)-1])
          #   if (class(t[2,2])=="character"){
          #     nr=dim(t)[[1]]
          #     nc=dim(t)[[2]]
          #     t=sapply(t, as.numeric)
          #     t=matrix(t,nr,nc)
          #   }
          #   incProgress(3/5, detail = paste("Step: Calculating"))
          #   newfit<-lmFit(t,design = modmat)
          #
          #   eb_newfit<-eBayes(newfit)
          #   incProgress(4/5, detail = paste("Step: Rendering table"))
          #   #names(sort(exp_genes_mean,decreasing = T)[1:input$nogenes])
          #   toptable_sample<-topTable(eb_newfit,number = ncol(alldat)-1)
          # }) ## end of calculating DE progress

          output$GeneOntology <- renderTable({
            rownames(gene_ids_all)<-gene_ids_all[,3]
            withProgress(message = 'calculating Gene Ontology',detail=NULL, value = 0, {
              incProgress(1/3, detail = paste("Step: Getting gene IDs"))
              degenes_ids<-degenes_table_ord[1:input$nogenes,1]
              degenes_int<-intersect(degenes_ids,rownames(gene_ids_all))
              degenes<-gene_ids_all[degenes_int,3]
              incProgress(2/3, detail = paste("Step: Calculating "))
              for (i in 1:nrow(go_dispdat)) {
                go_dispdat[i,2]<-length(Hs.c5[[i]])
                go_dispdat[i,3]<-length(intersect(degenes,Hs.c5[[i]]))
                go_dispdat[i,4]<-phyper(length(intersect(degenes,Hs.c5[[i]])),length(Hs.c5[[i]]),ncol(scdata_subset)-1-length(Hs.c5[[i]]),length(degenes),lower.tail = F)
                int_genes_ids<-intersect(degenes,Hs.c5[[i]])
                int_genes<-gene_ids_all[int_genes_ids,1]
                if (length(int_genes)>0){
                  go_dispdat[i,5]<-(paste(int_genes,collapse=", "))
                } else {
                  go_dispdat[i,5]<-as.character("0")
                }
              }
              incProgress(3/3, detail = paste("Step: Getting Geneontology"))
              go_ord<-go_dispdat[order(go_dispdat[,4]),]
              showNotification("GeneOntology calculation finished")
              output$downloadGO <- downloadHandler(
                filename = function() {
                  paste("GO_data", ".csv", sep = "")
                },
                content = function(file) {
                  write.csv(go_ord, file, row.names = FALSE)
                }
              )
              head(go_ord,n = 10)
            })
          },bordered = T)

          ### DE gene buttons implementation:
          DEgenes<-degenes_table_ord[1:input$nogenes,1] # a vector of characters
          for (i in 1:length(DEgenes))
          {
            assign("degenenames",c(degenenames,paste(DEgenes[i]," ",seq="")),envir=env)
          }
          output$DEbuttons <- renderUI({
            lapply(
              X = 1:length(DEgenes),
              FUN = function(i){
                actionButton(paste(DEgenes[i]," ",seq=""),paste(DEgenes[i]," ",seq=""))
              }
            )
          })
          ### maintain DE gene buttons

          lapply(
            X = 1:length(DEgenes),
            FUN = function(i){

              o<-observeEvent(input[[paste(DEgenes[i]," ",seq="")]], {
                showNotification(paste("showing ", DEgenes[i],"'s expression",sep=""),duration=5)
                output$plot <- renderPlotly({
                  plot_ly(
                    df,
                    #width=3000,height=4000,
                    x = df$x1, y = df$x2,
                    text = ~paste("label: ", as.factor(df$y)),
                    color = (scdata_subset[[DEgenes[i]]]),
                    key = row.names(df),
                    type = 'scatter',
                    mode = 'markers'
                  )%>% layout(dragmode = "lasso",title=paste("Value of ",DEgenes[i],sep=""))
                })
              }
              )

              assign("debuttons",c(debuttons,isolate(o)),envir =env)
            })
          ## end of maintaining buttons

          ##### constructing hgnc_filt using informations in the marker
          ENTREZID=array()
          SYMBOL=array()
          for (i in 1:length(markers)){
            ENTREZID=c(ENTREZID,markers[[as.character(i-1)]][["indices"]])
            SYMBOL=c(SYMBOL,markers[[as.character(i-1)]][["outp_names"]])
          }
          for (i in 1:length(SYMBOL)){
            if (is.na(SYMBOL[i])){
              SYMBOL[i]="NA"
            }
          }
          SYMBOL=data.frame(SYMBOL)
          ENTREZID=data.frame(ENTREZID)
          SYMBOL=distinct(SYMBOL)
          ENTREZID=distinct(ENTREZID)
          hgnc_filt=data.frame(SYMBOL,ENTREZID)
          row.names(hgnc_filt)=as.character(SYMBOL[[1]])
          ##### end of constructing hgnc_filt dataframe of genename,id

          ### KEGG panel
          output$KEGG <- renderTable({
            rownames(gene_ids_all)<-gene_ids_all[,3]
            withProgress(message = 'calculating KEGG',detail=NULL, value = 0, {
              incProgress(1/3, detail = paste("Step: Getting gene IDs"))
              degenes_ids<-degenes_table_ord[1:input$nogenes,1]
              degenes_int<-intersect(degenes_ids,rownames(gene_ids_all))
              degenes<-gene_ids_all[degenes_int,3]
              incProgress(2/3, detail = paste("Step: Calculating "))
              for (i in 1:nrow(kegg_dispdat)) {
                kegg_dispdat[i,2]<-length(kegg_genelists[[i]])
                kegg_dispdat[i,3]<-length(intersect(degenes,kegg_genelists[[i]]))
                kegg_dispdat[i,4]<-phyper(length(intersect(degenes,kegg_genelists[[i]])),length(kegg_genelists[[i]]),ncol(scdata_subset)-1-length(kegg_genelists[[i]]),length(degenes),lower.tail = F)
                int_genes_ids<-intersect(degenes,kegg_genelists[[i]])
                int_genes<-gene_ids_all[int_genes_ids,1]
                if (length(int_genes)>0){
                  kegg_dispdat[i,5]<-(paste(int_genes,collapse=", "))
                } else {
                  kegg_dispdat[i,5]<-as.character("0")
                }
              }
              incProgress(3/3, detail = paste("Step: Formating"))
              kegg_ord<-kegg_dispdat[order(kegg_dispdat[,4]),]
              showNotification("KEGG calculation finished")
              output$downloadKEGG <- downloadHandler(
                filename = function() {
                  paste("KEGG_data", ".csv", sep = "")
                },
                content = function(file) {
                  write.csv(kegg_ord, file, row.names = FALSE)
                }
              )
              head(kegg_ord,n=10)
            })
          },bordered = T)

          ### Markers panel
          output$Markers <- renderTable({
            rownames(gene_ids_all)<-gene_ids_all[,1]
            withProgress(message = 'calculating Markers Intersect',detail=NULL, value = 0, {
              incProgress(1/3, detail = paste("Step: Getting gene IDs"))
              degenes_new<-degenes_table_ord[1:input$nogenes,1]
              degenes<-intersect(degenes_new,rownames(gene_ids_all))
              #degenes<-gene_ids_all[degenes_int,1]
              incProgress(2/3, detail = paste("Step: Calculating hypergeom"))
              for (i in 1:nrow(hypergeom)) {
                hypergeom[i,1]<-names(markers_genelists_list)[i]
                hypergeom[i,2]<-length(markers_genelists_list[[i]])
                hypergeom[i,3]<-length(intersect(degenes,markers_genelists_list[[i]]))
                hypergeom[i,4]<-phyper(length(intersect(degenes,markers_genelists_list[[i]])),length(markers_genelists_list[[i]]),ncol(scdata_subset)-1-length(markers_genelists_list[[i]]),length(degenes),lower.tail = F)
                int_genes_ids<-intersect(degenes,markers_genelists_list[[i]])
                int_genes<-gene_ids_all[int_genes_ids,1]
                if (length(int_genes)>0){
                  hypergeom[i,5]<-(paste(int_genes,collapse=", "))
                } else {
                  hypergeom[i,5]<-as.character("0")
                }
              }
              incProgress(3/3, detail = paste("Step: Presenting"))
              hypergeom_ord<-hypergeom[order(hypergeom[,4]),]
              showNotification("Markers Intersect calculation finished")
              output$downloadMKS <- downloadHandler(
                filename = function() {
                  paste("Markers_data", ".csv", sep = "")
                },
                content = function(file) {
                  write.csv(hypergeom_ord, file, row.names = FALSE)
                }
              )
              hypergeom_ord[1:10,]
            })
          },bordered = T)

          #Top expressed genes
          #output$topgenes<-renderPrint({
          # aveeexpr<-colMeans(selecteddat[,2:ncol(selecteddat)],na.rm = T)
          #aveexpr_sort<-sort(aveexpr,decreasing = T)
          #aveexpr[1:input$nogenes]
          #})
          ### Msigdb panel
          output$Msigdb <- renderTable({
            rownames(gene_ids_all)<-gene_ids_all[,3]
            withProgress(message = 'calculating Msigdb',detail=NULL, value = 0, {
              incProgress(1/3, detail = paste("Step: Getting gene IDs"))
              degenes_ids<-degenes_table_ord[1:input$nogenes,1]
              degenes_int<-intersect(degenes_ids,rownames(gene_ids_all))
              degenes<-gene_ids_all[degenes_int,3]
              incProgress(2/3, detail = paste("Step: Calculating"))
              for (i in 1:nrow(msig_dispdat)) {
                msig_dispdat[i,2]<-length(Hs.c2[[i]])
                msig_dispdat[i,3]<-length(intersect(degenes,Hs.c2[[i]]))
                msig_dispdat[i,4]<-phyper(length(intersect(degenes,Hs.c2[[i]])),length(Hs.c2[[i]]),ncol(scdata_subset)-1-length(Hs.c2[[i]]),length(degenes),lower.tail = F)
                int_genes_ids<-intersect(degenes,Hs.c2[[i]])
                int_genes<-gene_ids_all[int_genes_ids,1]
                if (length(int_genes)>0){
                  msig_dispdat[i,5]<-(paste(int_genes,collapse=", "))
                } else {
                  msig_dispdat[i,5]<-as.character("0")
                }
              }
              incProgress(3/3, detail = paste("Step: Presenting"))
              msig_ord<-msig_dispdat[order(msig_dispdat[,4]),]
              showNotification("Msigdb calculation finished")
              output$downloadMSIG <- downloadHandler(
                filename = function() {
                  paste("MsigDB_data", ".csv", sep = "")
                },
                content = function(file) {
                  write.csv(msig_ord, file, row.names = FALSE)
                }
              )
              msig_ord[1:10,]

            })
          },bordered = T)
          degenes_table_ord[1:input$nogenes,]
        })
      })   ###end of get de gene button



      # selecteddat=NULL
      # flag=0
      # if (input$getdegenes>set){
      # if (input$getdegenes==0)
      # {
      #   return()
      # }
      # if (length(debuttons)>0){
      #   for (i in 1:length(debuttons)){
      #     debuttons[[i]]$destroy()
      #     #debuttons[[1]]<-NULL
      #   }
      # }
      # assign("debuttons",NULL,envir=env)  # disable previous buttons

      #assign("set", set+1, envir = env)


      #flag=1
    }

    # if (input$DEsubsets>sets){
    #     if (input$DEsubsets==0)
    #     {
    #         return()
    #     }




    #assign("sets", sets+1, envir = env)


    #     flag=1
    # }



    # cluster_selection: idxes<-which(expr_data$cluster==as.numeric(input$cluforanalysis))


    ### finisehd calculating selected data and rest data
    ### start calculating DE genes


    ############ No loNger used
    ###############################################################BOTTOM OF MAIN PANEL: (CLUSTERS & INTERSECTIONS)
    # Adding tabset panel corresponds to each Cluster
    # (at the bottom of the main panel)

    # for (i in 1:length(names(markers))) {
    #   if (i == 1) {
    #     insertTab(
    #       inputId = "switcher",
    #       position="after",
    #       tabPanel(
    #         as.character(i-1),
    #         paste("Cluster", as.character(i-1), " intersections", seq=""),
    #         tags$div(id = paste('placeholder', as.character(i-1), sep=""))
    #       ),
    #       target= "No selection"
    #     )
    #   } else {
    #     insertTab(
    #       inputId = "switcher",
    #       position = "after",
    #       tabPanel(
    #         as.character(i-1),
    #         paste("Cluster",as.character(i-1)," intersections",seq=""),
    #         tags$div(id = paste('placeholder',as.character(i-1),sep=""))),
    #       target = as.character(i-2)
    #     )
    #   }
    # }
    #
    # ### update tabset panel according to the point selected
    # o<-observeEvent(event_data("plotly_click", priority = "event"), {
    #   cldat <- event_data("plotly_click")
    #   selectedi = "No selection"
    #   selectedi <- filter(df, round(x1, 5) == round(cldat$x, 5))$y
    #   updateTabsetPanel(session, "switcher", selected=as.character(selectedi))
    #   i = df[, 3][cldat$pointNumber]
    #
    #   # if (previous_i!=-1){
    #   #   for (j in 1:length(c_intersections[[i]])){
    #   #     removeUI(
    #   #       selector = paste("#",c_intersections[[i]][j],"")
    #   #     )
    #   #   }
    #   # }
    # })
    #
    #
    #
    # assign("switcher",c(switcher,o),envir = env)
    # ##############  Adding intersection buttons to corresponding tab panel
    # # get the marker list intersections
    #
    # # step 1 c_intersection is a list of intersections in each cluster
    # c_intersections <- list("")
    # assign("clusters",length(names(markers)),envir = env)
    #
    # for (i in 1:length(names(markers))) {
    #   intersection <- markers[[as.character(i - 1)]][["lvl1_intersec"]]
    #   #intersection <- strsplit(str_replace_all(
    #   #           intersection,"([\\[])|([\\]])|([\n])|([\\'])", "")," ")
    #   if (isEmpty(intersection)){       ##intersection could be "numeric(0)" if there is no intersection in this cluster
    #     intersection=" "
    #   }
    #   intersection <- list(intersection)
    #   c_intersections <- append(c_intersections, intersection)
    # }
    # c_intersections[1] <- NULL
    #
    # #step 2 total_intersection is the total gene in the intersection
    # total_intersections = c()
    # c_seen = c()
    #
    # for (i in 1:length(c_intersections)) {
    #   for (j in 1:length(c_intersections[[i]])) {
    #     total_intersections <- c(
    #       total_intersections, c_intersections[[i]][j])
    #   }
    # }
    #
    # #step3 add a symbol to gene names that appear twice.
    # #so we won't have duplicated button IDs later
    # flag <- 0
    # tmp = c()
    # c_updated <- c()
    # for (x in 1:length(c_intersections)) {
    #   for (k in 1:length(c_intersections)) {
    #     for (i in 1:length(c_intersections[[k]])) {
    #       for (j in 1:length(total_intersections)) {
    #         if ((c_intersections[[k]][i] %in% c_seen) == FALSE ){
    #           c_seen <- c(c_seen, c_intersections[[k]][i])
    #           c_updated <- c(c_updated, k * 1000 + i)
    #         }
    #         if (identical(
    #           c_intersections[[k]][i], total_intersections[j])
    #           && (c_intersections[[k]][i] %in% c_seen)
    #           && (((k*1000+i) %in% c_updated)==FALSE)) {
    #           c_intersections[[k]][i] <- paste(
    #             c_intersections[[k]][i],"-",sep="")
    #           total_intersections <- c(
    #             total_intersections, c_intersections[[k]][i])
    #           c_seen <- c(c_seen, c_intersections[[k]][i])
    #           c_updated <- c(c_updated, k * 1000 + i)
    #         }
    #       }
    #     }
    #   }
    # }
    #
    # ##Step4: Adding buttons into corresponding tabpanels
    # for (i in 1:length(c_intersections)){
    #   for (j in 1:length(c_intersections[[i]])){
    #     textt<-c_intersections[[i]][j]
    #
    #     for (k in 1:length(total_intersections)){
    #       if (identical(
    #         total_intersections[k],
    #         as.character(strsplit(c_intersections[[i]][j], "-")))) {
    #         #showNotification(textt,duration=NULL)
    #         textt <- total_intersections[k]
    #         #showNotification(textt,duration=NULL)
    #         break
    #       }
    #     }
    #     if (substr(c_intersections[[i]][j],1,1)!=" "){
    #       insertUI(
    #         selector = paste("#placeholder",as.character(i-1),sep=""),
    #         #where = "afterEnd",
    #         ui = actionButton(c_intersections[[i]][j],textt)
    #       )
    #     }
    #   }
    # }
    #
    # ## Final step: Maintaining those buttons:
    # lapply(
    #   X = 1:length(total_intersections),
    #   FUN = function(i) {
    #     o<-observeEvent(input[[total_intersections[i]]], {
    #       rr = i
    #       for (j in 1:length(total_intersections)) {
    #         if (identical(total_intersections[j],
    #                       as.character(strsplit(total_intersections[i],"-")))) {
    #           rr <- j
    #           break
    #         }
    #       }
    #       showNotification(
    #         paste("showing ", total_intersections[rr],
    #               "'s expression",sep = ""),
    #         duration = 5
    #       )
    #       output$plot <- renderPlotly({
    #         plot_ly(
    #           df, x = df$x1, y = df$x2,
    #           text = ~paste("label: ", as.factor(df$y)),
    #           color = expr_data[[total_intersections[rr]]],
    #           key = row.names(df),
    #           type = 'scatter',
    #           mode = 'markers'
    #         ) %>% layout(dragmode = "lasso",
    #                      title = paste(
    #                        "Value of ", total_intersections[rr], sep=""))
    #       })
    #       observeEvent(input$showcard, {
    #         output$inc <- renderUI({
    #           x <- input$test
    #           getPage(total_intersections[rr])
    #         })
    #       })
    #     })
    #     assign("intersectbuttons",c(intersectbuttons,o),envir = env)
    #   })

    #### cluster analysis





    ###end of intersections



    ########################################START OF CHANGING CLUSTERS NAMES
    # updateSelectInput(session = session,
    #                   inputId = "chgcluster",
    #                   label = "Cluster to Rename",
    #                   choices = c(1:(length(c_intersections)-1)),
    #                   #selected = NULL)
    # )
    # observeEvent(input$chg,{
    #   removeUI(selector= paste("#placeholder",as.character(as.numeric(input$chgcluster)-1),sep=""))
    #
    #   insertTab(
    #     inputId = "switcher",
    #     position = "after",
    #     tabPanel(
    #       input$newcluster,
    #       paste("Cluster ",input$newcluster," intersections",seq=""),
    #       tags$div(id = paste("#placeholder",as.character(as.numeric(input$chgcluster)-1),sep=""))),
    #     target = as.character(length(c_intersections)-2)
    #   )
    #   i=length(c_intersections)
    #   for (j in 1:length(c_intersections[[i]])){
    #     textt<-c_intersections[[i]][j]
    #
    #     for (k in 1:length(total_intersections)){
    #       if (identical(
    #         total_intersections[k],
    #         as.character(strsplit(c_intersections[[i]][j], "-")))) {
    #         #showNotification(textt,duration=NULL)
    #         textt <- total_intersections[k]
    #         #showNotification(textt,duration=NULL)
    #         break
    #       }
    #     }
    #     if (substr(c_intersections[[i]][j],1,1)!=" "){
    #       insertUI(
    #         selector = paste("#placeholder",as.character(i-1),sep=""),
    #         #where = "afterEnd",
    #         ui = actionButton(c_intersections[[i]][j],textt)
    #       )
    #     }
    #   }
    # })

    #########################################END OF CHANGING CLUSTERS NAMES












  })




  ####################################    END OF RUN CURRENT CONFIG


})


#   #showNotification(as.character(round(df$x1[[1]],5)),duration=NULL)
# #UPDATE
# ############################################################################
# observeEvent(input$update, {
#
#   plotid<-plotid+1
#  # withReactiveDomain(as.character(plotid),{
#
#   df<-NULL
#   for (i in length(total_intersections)){
#     removeUI(selector=paste("#",total_intersections,sep=""))
#   }
#   for (i in 1:clusters){
#     removeTab(inputId="switcher", target=as.character(i-1))
#   }
#
#   # insertUI(selector="#placeholder",where="afterEnd",
#   #          tabsetPanel(
#   #            id = "switcher",
#   #            tabPanel("No selection", "No selection"))
#   # )
#   df <- isolate(get_plot_data())
#   #get the gene expression data
#   expr_data=matrix(pipe$x,ncol=length(pipe$col_ids),dimnames=list(1:length(pipe$x[,1]),pipe$col_ids))
#   expr_data=data.frame(expr_data)
#   expr_data$cluster=df$y
#   updateSelectInput(session=session, inputId="color", label = "Select colour value:", choices = c("cluster",names(expr_data)),
#                     selected = NULL)
#   updateSelectInput(session=session, inputId="newlabels", label = "Select label", choices =levels(as.factor(expr_data[,length(expr_data)])),selected=NULL)
#   markers=pipe$markers
#   #Adding tabset panel corresponds to each cluster
#   for (i in 1:length(names(markers))){
#     if (i==1)
#     {
#       insertTab(inputId = "switcher", position="after",
#                 tabPanel(as.character(i-1), paste("Cluster",as.character(i-1)," intersections",seq=""),tags$div(id = paste('placeholder',as.character(i-1),sep=""))),
#                 target= "No selection"
#       )
#     }
#     else
#     {
#       insertTab(inputId = "switcher", position="after",
#                 tabPanel(as.character(i-1), paste("Cluster",as.character(i-1)," intersections",seq=""),tags$div(id = paste('placeholder',as.character(i-1),sep=""))),
#                 target= as.character(i-2)
#       )
#     }
#   }
#
#
#   ##Adding intersection buttons
#   # step 1 c_intersection is a list of intersections in each cluster
#   c_intersections <- list("")
#   clusters<-length(names(markers))
#   for (i in 1:length(names(markers))){
#     intersection <- markers[[as.character(i-1)]][["lvl1_intersec"]]
#     #intersection <- strsplit(str_replace_all(intersection,"([\\[])|([\\]])|([\n])|([\\'])", "")," ")
#     intersection <- list(intersection)
#     c_intersections <- append(c_intersections,intersection)
#   }
#   c_intersections[1]<-NULL
#   #step 2 total_intersection is the total gene in the intersection
#   total_intersections=c()
#   c_seen=c()
#   for (i in 1:length(c_intersections)){
#     for (j in 1:length(c_intersections[[i]])){
#       total_intersections <- c(total_intersections,c_intersections[[i]][j])
#     }
#   }
#   #step3 add a symbol to gene names that appear twice.
#   #so we won't have duplicated button IDs later
#   flag<-0
#   tmp=c()
#   c_updated <- c()
#   for (x in 1:length(c_intersections)){
#     for (k in 1:length(c_intersections)){
#       for (i in 1:length(c_intersections[[k]])){
#         for (j in 1:length(total_intersections)){
#           if ((c_intersections[[k]][i] %in% c_seen)==FALSE){
#             c_seen<-c(c_seen,c_intersections[[k]][i])
#             c_updated<-c(c_updated,k*1000+i)
#           }
#           if (identical(c_intersections[[k]][i] , total_intersections[j]) && (c_intersections[[k]][i] %in% c_seen) && (((k*1000+i) %in% c_updated)==FALSE)){
#             c_intersections[[k]][i]<-paste(c_intersections[[k]][i],"-",sep="")
#             total_intersections<-c(total_intersections,c_intersections[[k]][i])
#             c_seen<-c(c_seen,c_intersections[[k]][i])
#             c_updated<-c(c_updated,k*1000+i)
#           }
#         }
#       }
#     }
#   }
#   ##Step4: Adding buttons into corresponding tabpanels
#   for (i in 1:length(c_intersections)){
#     for (j in 1:length(c_intersections[[i]])){
#       textt<-c_intersections[[i]][j]
#
#       for (k in 1:length(total_intersections)){
#         if (identical(total_intersections[k],as.character(strsplit(c_intersections[[i]][j],"-")))){
#           #showNotification(textt,duration=NULL)
#           textt<-total_intersections[k]
#           #showNotification(textt,duration=NULL)
#           break
#         }
#       }
#       insertUI(
#         selector = paste("#placeholder",as.character(i-1),sep=""),
#         #where = "afterEnd",
#         ui = actionButton(c_intersections[[i]][j],textt)
#       )
#     }
#   }
#   ## Final step: Maintaining those buttons:
#   lapply(
#     X = 1:length(total_intersections),
#     FUN = function(i){
#       observeEvent(input[[total_intersections[i]]], {
#
#         rr=i
#         for (j in 1:length(total_intersections)){
#           if (identical(total_intersections[j],as.character(strsplit(total_intersections[i],"-")))){
#             rr<-j
#             break
#           }
#         }
#
#         #showNotification(paste("showing ", total_intersections[rr],"'s expression",sep=""),duration=5)
#         output$plot <- renderPlotly({
#           plot_ly(
#             df,
#             x = df$x1, y = df$x2,
#             text = ~paste("label: ", as.factor(df$y)),
#             #color = as.factor(df$y)
#             color = expr_data[[total_intersections[rr]]],
#             key = row.names(df)
#           ) %>% layout(dragmode = "lasso")
#         })
#         observeEvent(input$showcard,{
#           output$inc<-renderUI({
#             x <- input$test
#             getPage(total_intersections[rr])
#           })
#           #showNotification(paste("showing ", total_intersections[rr],"'s gene card",sep=""),duration=5)
#         })
#       })
#     }
#   )
#
#
#   updateSelectInput(session=session, inputId="color", label= "Select color value" ,choices = names(expr_data),
#                     selected = NULL)                                                                   ### and update color values that can be selected
#
#   output$plot <- renderPlotly({
#     plot_ly(
#       df,
#       x = df$x1, y = df$x2,
#       text = ~paste("label: ", as.factor(df$y)),
#       color = as.factor(df$y)
#     ) %>% layout(dragmode = "lasso")
#   })
#   observeEvent(event_data("plotly_click",priority = "event"), {
#     cldat <- event_data("plotly_click")
#
#     #clnewdat<-sc_data_newclusts$labels
#     #plotdata[,3]
#     selectedi="No selection"
#     selectedi <- filter(df,round(x1,5)==round(cldat$x,5))$y
#     showNotification(as.character(round(cldat$x,5)),duration=NULL)
#     #showNotification(as.character(round(df$x1[[1]],5)),duration=NULL)
#     #showNotification(as.character(selectedi),duration=NULL)
#     updateTabsetPanel(session, "switcher", selected = as.character(selectedi))
#
#     i=df[,3][cldat$pointNumber]
#   })
#   ##DE GENE IMPLEMENTATION
#   scdata_subset=expr_data
#   observeEvent(input$getdegenes,{
#     output$genes <- renderPrint({
#       d <- event_data("plotly_selected")
#       selecteddat<-scdata_subset[as.numeric(d$key),2:ncol(scdata_subset)]
#       restdat<-scdata_subset[-as.numeric(d$key),2:ncol(scdata_subset)]
#       #exp_genes_mean<-colSums(exp_genes)/nrow(exp_genes)
#       labelsdat<-as.factor(c(rep("selected",nrow(selecteddat)),rep("notselected",nrow(restdat))))
#       alldat<-rbind(selecteddat,restdat)
#       alldat<-data.frame(alldat,labelsdat)
#       modmat<-model.matrix(~labelsdat,data = alldat)
#       newfit<-lmFit(t(alldat[,1:ncol(alldat)-1]),design = modmat)
#       eb_newfit<-eBayes(newfit)
#       #names(sort(exp_genes_mean,decreasing = T)[1:input$nogenes])
#       toptable_sample<-topTable(eb_newfit,number = ncol(alldat)-1)
#       output$GeneOntology <- renderPrint({
#         geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
#         gotable<-goana(geneids)
#         gotable[1:10,]
#       })
#       # DE gene buttons implementation:
#       DEgenes<-rownames(toptable_sample[1:input$nogenes,]) # a vector of characters
#       output$DEbuttons <- renderUI({
#         lapply(
#           X = 1:length(DEgenes),
#           FUN = function(i){
#             actionButton(paste(DEgenes[i]," ",seq=""),paste(DEgenes[i]," ",seq=""))
#           }
#         )
#       })
#       lapply(
#         X = 1:length(DEgenes),
#         FUN = function(i){
#           observeEvent(input[[paste(DEgenes[i]," ",seq="")]], {
#             showNotification(paste("showing ", DEgenes[i],"'s expression",sep=""),duration=5)
#             output$plot <- renderPlotly({
#               plot_ly(
#                 df,
#                 x = df$x1, y = df$x2,
#                 text = ~paste("label: ", as.factor(df$y)),
#                 color = (scdata_subset[[DEgenes[i]]]),
#                 key = row.names(df)
#               )%>% layout(dragmode = "lasso")
#             })
#           }
#           )
#         })
#       output$KEGG <- renderPrint({
#         geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
#         keggtable<-kegga(geneids)
#         keggtable[1:10,]
#       })
#       output$Markers <- renderPrint({
#         degenes<-rownames(toptable_sample[1:input$nogenes,])
#         hypergeom[1,2]<-phyper(length(intersect(degenes,c_intersections[[1]])),length(c_intersections[[1]]),ncol(scdata_subset)-1-length(c_intersections[[1]]),length(degenes),lower.tail = F)
#         hypergeom[2,2]<-phyper(length(intersect(degenes,c_intersections[[2]])),length(c_intersections[[2]]),ncol(scdata_subset)-1-length(c_intersections[[2]]),length(degenes),lower.tail = F)
#         hypergeom[3,2]<-phyper(length(intersect(degenes,c_intersections[[3]])),length(c_intersections[[3]]),ncol(scdata_subset)-1-length(c_intersections[[3]]),length(degenes),lower.tail = F)
#         hypergeom[4,2]<-phyper(length(intersect(degenes,c_intersections[[4]])),length(c_intersections[[4]]),ncol(scdata_subset)-1-length(c_intersections[[4]]),length(degenes),lower.tail = F)
#         hypergeom
#       })
#       toptable_sample[1:input$nogenes,]
#     })
#   })
# ############################################################################
# }#end of obs event expression
# )# end of observe event
