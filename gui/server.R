library(reticulate)
library(shiny)
library(plotly)
library(ggplot2)
library(stringr)
library(limma)
library(GO.db)

#ids <- read.csv('datasets/spleen/spleen.csv', nrows = 1, header = FALSE)
#X <- read.csv('datasets/spleen/spleen.csv', skip = 1, header = FALSE)

# Load python

#use_virtualenv("modAL")
#use_virtualenv("py3.7")
#use_python("D")
source_python('__init__.py')

# Create pipeline
pipe <- Pipeline()

## getpage function for genecard
getPage<-function(genename) {
  url <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",genename,sep="")
  return(browseURL(url))
}

intersect<-function (x, y) 
{
  y <- as.vector(y)
  unique(y[match(as.vector(x), y, 0L)])
}


server <- shinyServer(function(input, output, session) {
  #catch file
  observeEvent(input$file1,{
  req(input$file1)
  tryCatch(
    {
      f <- read.csv(input$file1$datapath,
                     #header = input$header,
                     #sep = input$sep,
                     #quote = input$quote)
          )
      fname=strsplit(input$file1$name,".",fixed=TRUE)[[1]][1]
      dir.create(paste(getwd(),"/datasets/",fname,sep=""))
      write.csv(f,paste(getwd(),"/datasets/",fname,"/",fname,sep=""))
      showNotification("File uploaded. Now you can specify the new dataset in the configuration")
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    }
    
    #js$reset()
  )
    
  })
  
  
  
  #start of data fetching and processing
  observeEvent(input$reset, {js$reset()})
  plotid=0
  #running clustering
  ############################################################################
  get_plot_data <- isolate(function(){
    dim_method <- input$dim_method
    dim_n_components <- input$dim_n_components
    clu_method <- input$clu_method
    eval_method <- input$eval_method
    clu_n_clusters <- input$clu_n_clusters
    clu_n_jobs <- input$clu_n_jobs
    mark_method <- 'TTest'
    mark_alpha <- input$mark_alpha
    mark_markers_n <- input$mark_markers_n
    mark_correction <- input$mark_correction
    mark_n_jobs <- input$mark_n_jobs
    con_method <- 'Converter'
    con_convention <- input$con_convention
    #con_path <- input$con_path
    con_path <- 'markers/gene_id_name.csv'
    ide_method <- 'HyperGeom'
    #ide_path <- input$ide_path
    ide_path <- 'markers/cell_type_marker.json'
    ide_tissue <- input$ide_tissue
    vis_method <- input$vis_method
    ssc_method <- input$ssc_method

    # pipe$run(dim_method=dim_method, dim_n_components=dim_n_components,
    #         clu_method=clu_method, eval_method=eval_method,
    #         clu_n_clusters=clu_n_clusters, ssc_method=ssc_method)

    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 6

      incProgress(1/n, detail = paste("Step: PCA"))
      x_emb = pipe$get_emb(
            method=dim_method, n_components=dim_n_components)

      incProgress(1/n, detail = paste("Step: Clustering"))
      labels <- pipe$get_labels(
            method=clu_method,
            eval_method=eval_method, n_clusters=clu_n_clusters,
            n_jobs=clu_n_jobs)

      incProgress(1/n, detail = paste("Step: Finding markers"))
      markers <- pipe$get_markers(
            method='TTest',
            alpha=mark_alpha, markers_n=mark_markers_n,
            correction=mark_correction, n_jobs=mark_n_jobs)

      incProgress(1/n, detail = paste("Step: Converting names"))
      markers <- pipe$convert(
            method=con_method,
            convention=con_convention, path=con_path)

      incProgress(1/n, detail = paste("Step: Identifying cells"))
      markers <- pipe$identify(
            method="HyperGeom",
            path=ide_path, tissue=ide_tissue)

      incProgress(1/n, detail = paste("Step: Visualizing"))
      x_emb_2d <- pipe$get_emb_2d(x_emb, y=labels, method=vis_method)
      df <- data.frame(x1 = x_emb_2d[, 1],
                      x2 = x_emb_2d[, 2],
                      y = labels)
    })
    return(df)
  })
  df <- isolate(get_plot_data())
  markers<-pipe$markers
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
  
  #get the gene expression data
  expr_data=matrix(pipe$x,ncol=length(pipe$col_ids),dimnames=list(1:length(pipe$x[,1]),pipe$col_ids))
  expr_data=data.frame(expr_data)
  expr_data$cluster=df$y
  updateSelectInput(session=session, inputId="color", label = "Select colour value:", choices = c("cluster",names(expr_data)),
                    selected = NULL)
  updateSelectInput(session=session, inputId="newlabels", label = "Select label", choices =levels(as.factor(expr_data[,length(expr_data)])),selected=NULL)
  markers=pipe$markers
  
  #create object to store hypergeometric marker results
  markers<-c("Blood - CD1C+ B dendritic cell","Kidney - Cancer Stem cell","Liver - CD4+ cytotoxic T cell","Kidney - ErythroBlast")
  pvals<-double(length = 4)
  hypergeom<-data.frame(markers,pvals)
  
  #Adding tabset panel corresponds to each cluster
  
  for (i in 1:length(names(markers))){
    if (i==1)
    {
        insertTab(inputId = "switcher", position="after",
              tabPanel(as.character(i-1), paste("Cluster",as.character(i-1)," intersections",seq=""),tags$div(id = paste('placeholder',as.character(i-1),sep=""))),
              target= "No selection"
              )
    }
    else
    {
      insertTab(inputId = "switcher", position="after",
                tabPanel(as.character(i-1), paste("Cluster",as.character(i-1)," intersections",seq=""),tags$div(id = paste('placeholder',as.character(i-1),sep=""))),
                target= as.character(i-2)
      )
    }
  } 
  
 
  ##Adding intersection buttons 
  # step 1 c_intersection is a list of intersections in each cluster
  c_intersections <- list("")
  clusters<-length(names(markers))
  for (i in 1:length(names(markers))){
    intersection <- markers[[as.character(i-1)]][["lvl1_intersec"]]
    #intersection <- strsplit(str_replace_all(intersection,"([\\[])|([\\]])|([\n])|([\\'])", "")," ")
    intersection <- list(intersection)
    c_intersections <- append(c_intersections,intersection)
  }
  c_intersections[1]<-NULL
  #step 2 total_intersection is the total gene in the intersection
  total_intersections=c()
  c_seen=c()
  for (i in 1:length(c_intersections)){
    for (j in 1:length(c_intersections[[i]])){
      total_intersections <- c(total_intersections,c_intersections[[i]][j])
    }
  }
  #step3 add a symbol to gene names that appear twice.
  #so we won't have duplicated button IDs later
  flag<-0
  tmp=c()
  c_updated <- c()
  for (x in 1:length(c_intersections)){
    for (k in 1:length(c_intersections)){
      for (i in 1:length(c_intersections[[k]])){
        for (j in 1:length(total_intersections)){
          if ((c_intersections[[k]][i] %in% c_seen)==FALSE){
            c_seen<-c(c_seen,c_intersections[[k]][i])
            c_updated<-c(c_updated,k*1000+i)
          }
          if (identical(c_intersections[[k]][i] , total_intersections[j]) && (c_intersections[[k]][i] %in% c_seen) && (((k*1000+i) %in% c_updated)==FALSE)){
            c_intersections[[k]][i]<-paste(c_intersections[[k]][i],"-",sep="")
            total_intersections<-c(total_intersections,c_intersections[[k]][i])
            c_seen<-c(c_seen,c_intersections[[k]][i])
            c_updated<-c(c_updated,k*1000+i)
          }
        }
      }
    }
  }
  ##Step4: Adding buttons into corresponding tabpanels 
  for (i in 1:length(c_intersections)){
    for (j in 1:length(c_intersections[[i]])){
      textt<-c_intersections[[i]][j]
      
      for (k in 1:length(total_intersections)){
        if (identical(total_intersections[k],as.character(strsplit(c_intersections[[i]][j],"-")))){
          #showNotification(textt,duration=NULL)
          textt<-total_intersections[k]
          #showNotification(textt,duration=NULL)
          break
        }
      }
      insertUI(
        selector = paste("#placeholder",as.character(i-1),sep=""),
        #where = "afterEnd",
        ui = actionButton(c_intersections[[i]][j],textt)
      )
    }
  }
  ## Final step: Maintaining those buttons:
  lapply(
    X = 1:length(total_intersections),
    FUN = function(i){
      observeEvent(input[[total_intersections[i]]], {
        
        rr=i
        for (j in 1:length(total_intersections)){
          if (identical(total_intersections[j],as.character(strsplit(total_intersections[i],"-")))){
            rr<-j
            break
          }
        }
        
        showNotification(paste("showing ", total_intersections[rr],"'s expression",sep=""),duration=5)
        output$plot <- renderPlotly({
          plot_ly(
            df,
            x = df$x1, y = df$x2,
            text = ~paste("label: ", as.factor(df$y)),
            #color = as.factor(df$y)
            color = expr_data[[total_intersections[rr]]],
            key = row.names(df)
          ) %>% layout(dragmode = "lasso")
        })
        observeEvent(input$showcard,{
          output$inc<-renderUI({
            x <- input$test
            getPage(total_intersections[rr])
          })
          #showNotification(paste("showing ", total_intersections[rr],"'s gene card",sep=""),duration=5)
        })
      })
    }
  )
  
  
  ###################################################################
  #gene card
  observeEvent(input$search, {
    if (input$searchgene %in% names(expr_data)){
      getPage(input$searchgene)
    }
    else{
      showNotification("Gene name does not exit")
    }
  })
  #gene card
  
  
  #run default plot
  ##############################################################
  output$plot <- renderPlotly({
    
    #factorize cluster labels (discrete instead of continuous)
    if (input$color == "cluster"){
      plotcols = as.factor(expr_data[[input$color]])
      
    } else {
      plotcols = expr_data[[input$color]]
      
    }
    
    plot_ly(
      df,
      x = df$x1, y = df$x2,
      text = ~paste("label: ", as.factor(df$y)),
      #color = as.factor(df$y)
      color = plotcols,
      key = row.names(df)
    ) %>% layout(dragmode = "lasso")
  })
  ############################################################################
  
  newlabs<-df[,3]
  names(newlabs)<-rownames(df)
  #plotdata<-data.frame(plotdata)
  #plotdata$key<-row.names(plotdata)
  #labeldats<-levels(as.factor(pipe$labels))
  labeldats<-levels(as.factor(df[,3]))
  
  observe({
    x <- input$text
    #Can use character(0) to remove all choices
    #labeldats<-reactive({value = labchoices})
    observeEvent(input$labeladd, {
      labeldats<<-union(labeldats,x)
      # Can also set the label and select items
      updateSelectInput(session, "newlabels",
                        label = paste("Select input label", length(x)),
                        choices = labeldats
      )
    })
  })
  
  #Title of the plot
  Title <- reactive({
    paste("Value of ", input$gene )
  })
  
  # output caption to ptibt title
  output$caption <- renderText({
    Title()
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
  
    observeEvent(input$labelupd, {
      newlabs[d$key]<-as.integer(input$newlabels)
      
      #(newlabs[d$key])
      
      output$Plot2<-renderPlotly({
        
        plot_ly(
          df,
          x=df[,1],y=df[,2],
          text = ~paste("label: ", as.factor(newlabs)),
          color = as.factor(newlabs)
        )%>% layout(dragmode = "lasso")
        
      })
    })
  })
  output$cc <- reactive({
    event_data("plotly_click")
  })
    observeEvent(event_data("plotly_click",priority = "event"), {
      cldat <- event_data("plotly_click")
      
      #clnewdat<-sc_data_newclusts$labels
      #plotdata[,3]
      selectedi="No selection"
      selectedi <- filter(df,round(x1,5)==round(cldat$x,5))$y
      #showNotification(as.character(round(cldat$x,5)),duration=NULL)
      #showNotification(as.character(selectedi),duration=NULL)
      updateTabsetPanel(session, "switcher", selected = as.character(selectedi))
      
      i=df[,3][cldat$pointNumber]
      # if (previous_i!=-1){
      #   for (j in 1:length(c_intersections[[i]])){
      #     removeUI(
      #       selector = paste("#",c_intersections[[i]][j],"")
      #     )
      #   }
      # }
      #showNotification(as.character(clnewdat[cldat$pointNumber+1]), duration = NULL)
    })
    
    
    ##DE GENE IMPLEMENTATION
    scdata_subset=expr_data
    observeEvent(input$getdegenes,{
      output$genes <- renderPrint({
        d <- event_data("plotly_selected")
        selecteddat<-scdata_subset[as.numeric(d$key),2:ncol(scdata_subset)]
        restdat<-scdata_subset[-as.numeric(d$key),2:ncol(scdata_subset)]
        #exp_genes_mean<-colSums(exp_genes)/nrow(exp_genes)
        labelsdat<-as.factor(c(rep("selected",nrow(selecteddat)),rep("notselected",nrow(restdat))))
        alldat<-rbind(selecteddat,restdat)
        alldat<-data.frame(alldat,labelsdat)
        modmat<-model.matrix(~labelsdat,data = alldat)
        newfit<-lmFit(t(alldat[,1:ncol(alldat)-1]),design = modmat)
        eb_newfit<-eBayes(newfit)
        #names(sort(exp_genes_mean,decreasing = T)[1:input$nogenes])
        toptable_sample<-topTable(eb_newfit,number = ncol(alldat)-1)
        output$GeneOntology <- renderPrint({
          geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
          gotable<-goana(geneids)
          go_ord<-gotable[order(gotable$P.DE),]
          go_ord[1:10,]
        })
        # DE gene buttons implementation:
        DEgenes<-rownames(toptable_sample[1:input$nogenes,]) # a vector of characters
        output$DEbuttons <- renderUI({
          lapply(
            X = 1:length(DEgenes),
            FUN = function(i){
              actionButton(paste(DEgenes[i]," ",seq=""),paste(DEgenes[i]," ",seq=""))
            }
          )
        })
        lapply(
          X = 1:length(DEgenes),
          FUN = function(i){
            observeEvent(input[[paste(DEgenes[i]," ",seq="")]], {
              showNotification(paste("showing ", DEgenes[i],"'s expression",sep=""),duration=5)
              output$plot <- renderPlotly({
                plot_ly(
                  df,
                  x = df$x1, y = df$x2,
                  text = ~paste("label: ", as.factor(df$y)),
                  color = (scdata_subset[[DEgenes[i]]]),
                  key = row.names(df)
                )%>% layout(dragmode = "lasso")
              })
            }
            )
          })
        output$KEGG <- renderPrint({
          geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
          keggtable<-kegga(geneids)
          kegg_ord<-keggtable[order(keggtable$P.DE),]
          kegg_ord[1:10,]
        })
        output$Markers <- renderPrint({
          degenes<-rownames(toptable_sample[1:input$nogenes,])
          hypergeom[1,2]<-phyper(length(intersect(degenes,c_intersections[[1]])),length(c_intersections[[1]]),ncol(scdata_subset)-1-length(c_intersections[[1]]),length(degenes),lower.tail = F)
          hypergeom[2,2]<-phyper(length(intersect(degenes,c_intersections[[2]])),length(c_intersections[[2]]),ncol(scdata_subset)-1-length(c_intersections[[2]]),length(degenes),lower.tail = F)
          hypergeom[3,2]<-phyper(length(intersect(degenes,c_intersections[[3]])),length(c_intersections[[3]]),ncol(scdata_subset)-1-length(c_intersections[[3]]),length(degenes),lower.tail = F)
          hypergeom[4,2]<-phyper(length(intersect(degenes,c_intersections[[4]])),length(c_intersections[[4]]),ncol(scdata_subset)-1-length(c_intersections[[4]]),length(degenes),lower.tail = F)
          hypergeom_ord<-hypergeom[order(hypergeom$pvals),]
          hypergeom_ord
        })
        output$Msigdb <- renderPrint({
        degenes<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
        for (i in 1:nrow(msig_dispdat)) {
          msig_dispdat[i,2]<-phyper(length(intersect(degenes,Hs.c2[[i]])),length(Hs.c2[[i]]),ncol(scdata_subset)-1-length(Hs.c2[[i]]),length(degenes),lower.tail = F)
        }
        msig_ord<-msig_dispdat[order(msig_dispdat$msigdb_pvals),]
        msig_ord[1:10,]
        })
        toptable_sample[1:input$nogenes,]
      })
    })
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
