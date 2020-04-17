library(reticulate)
library(shiny)
library(plotly)
library(ggplot2)
library(stringr)
library(limma)
library(GO.db)
library(rjson)


source_python('__init__.py')
pipe <- Pipeline(x = "default")

load('gui/Hs.c2')

source("gui/functions.R")
source("gui/pipe_functions.R")

################################################################# server
server <- shinyServer(function(input, output, session) {

    # Toggling of menus
    shinyjs::onclick("togglemain", {
        shinyjs::toggle(id = "mainpanel", anim = TRUE);
        shinyjs::hide(id = "configuration", anim = TRUE)
    })

    shinyjs::onclick("toggleconfig", {
        shinyjs::toggle(id = "configuration", anim = TRUE);
        shinyjs::hide(id = "mainpanel", anim = TRUE)
    })

    # rerun the app if "run with new configuration button" pressed
    # (this can avoid observing previous events)
    observeEvent(input$runconfigbtn, {
        js$reset()
    })

    # Upload dataset
    observeEvent(input$file1, {
        req(input$file1)
        tryCatch({
            writeDataset(input$file1$datapath, input$file1$name)
        }, error = function(e) {
            stop(safeError(e))
        })
    })

    df <- isolate(runPipe(pipe, input))
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

    #Requires HS.c2 to be loaded.
    msigdb_categories <- names(Hs.c2)
    msigdb_pvals <- double(length = length(msigdb_categories))
    msig_dispdat <- data.frame(msigdb_categories, msigdb_pvals)
    newlabs <- df[, 3]
    names(newlabs) <- rownames(df)
    labeldats <- levels(as.factor(df[,3]))

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
                      choices = levels(as.factor(expr_data[,length(expr_data)])),
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


    # gene card
    observeEvent(input$search, {
        if (input$searchgene %in% names(expr_data)) {
            getPage(input$searchgene) # function in functions.R
        }
        else {
            showNotification("Gene name does not exist.")
        }
    })

    ############################################## MAIN PANEL

    ### run default plot
    output$plot <- renderPlotly({
        #factorize cluster labels (discrete instead of continuous)
        if (input$color == "cluster") {
            plotcols = as.factor(expr_data[[input$color]])
        } else {
            plotcols = expr_data[[input$color]]
        }
        plot_ly(
            df,
            x = df$x1,
            y = df$x2,
            text = ~paste("label: ", as.factor(df$y)),
            color = plotcols,
            key = row.names(df)
        ) %>% layout(dragmode = "lasso",
                     title = paste("Value of ", input$color, sep=""))
    })

    ### updated plot
    output$brush <- renderPrint({
        d <- event_data("plotly_selected")
        observeEvent(input$labelupd, {
            newlabs[d$key] <- as.integer(input$newlabels)
            output$Plot2 <- renderPlotly({
                plot_ly(
                    df,
                    x = df[, 1],
                    y = df[, 2],
                    text = ~paste("label: ", as.factor(newlabs)),
                    color = as.factor(newlabs)
                ) %>% layout(dragmode = "lasso",
                            title = paste("Value of ", input$labelupd, sep=""))
            })
        })
    })

  ############################################## DE GENE IMPLEMENTATION
  scdata_subset=expr_data
  assign("s1", NULL, envir = .GlobalEnv)
  assign("s2", NULL, envir = .GlobalEnv)
  observeEvent(input$subset1,{
    d <- event_data("plotly_selected")
    assign("s1", d, envir = .GlobalEnv)
    showNotification("subset1 stored")
  })
  observeEvent(input$subset2,{
    d <- event_data("plotly_selected")
    assign("s2", d, envir = .GlobalEnv)
    showNotification("subset2 stored")
  })
    assign("sets", 0, envir = .GlobalEnv)
    assign("set", 0, envir = .GlobalEnv)
    toListen <- reactive({
      list(input$getdegenes,input$DEsubsets)

    })

    observeEvent(toListen(),{
      selecteddat=NULL
      if (input$getdegenes>set){
        assign("set", sets+1, envir = .GlobalEnv)
        d <- event_data("plotly_selected")
        selecteddat<-scdata_subset[as.numeric(d$key),2:ncol(scdata_subset)]
        restdat<-scdata_subset[-as.numeric(d$key),2:ncol(scdata_subset)]
        showNotification("DE genes one subset")
      }
      if (input$DEsubsets>sets){
        assign("sets", set+1, envir = .GlobalEnv)

        selecteddat<-scdata_subset[as.numeric(s1$key),2:ncol(scdata_subset)]
        restdat<-scdata_subset[as.numeric(s2$key),2:ncol(scdata_subset)]
        showNotification("DE genes two subsets")
      }
      else
      {

      }
      output$genes <- renderPrint({
        withProgress(message = 'calculating DE genes', value = 0, {

          #exp_genes_mean<-colSums(exp_genes)/nrow(exp_genes)
          labelsdat<-as.factor(c(rep("selected",nrow(selecteddat)),rep("notselected",nrow(restdat))))
          alldat<-rbind(selecteddat,restdat)
          alldat<-data.frame(alldat,labelsdat)
          modmat<-model.matrix(~labelsdat,data = alldat)
          newfit<-lmFit(t(alldat[,1:ncol(alldat)-1]),design = modmat)
          eb_newfit<-eBayes(newfit)
          #names(sort(exp_genes_mean,decreasing = T)[1:input$nogenes])
          toptable_sample<-topTable(eb_newfit,number = ncol(alldat)-1)
        })

        output$GeneOntology <- renderPrint({
          geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
          gotable<-goana(geneids)
          go_ord<-gotable[order(gotable$P.DE),]
          go_ord[1:10,]
        })

        ### DE gene buttons implementation:
        output$deinfo <- renderUI({
          h3("DE gene information:")
        })
        DEgenes<-rownames(toptable_sample[1:input$nogenes,]) # a vector of characters
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
            observeEvent(input[[paste(DEgenes[i]," ",seq="")]], {
              showNotification(paste("showing ", DEgenes[i],"'s expression",sep=""),duration=5)
              output$plot <- renderPlotly({
                plot_ly(
                  df,
                  x = df$x1, y = df$x2,
                  text = ~paste("label: ", as.factor(df$y)),
                  color = (scdata_subset[[DEgenes[i]]]),
                  key = row.names(df)
                )%>% layout(dragmode = "lasso",title=paste("Value of ",DEgenes[i],sep=""))
              })
            }
            )
          })
      ## end of maintaining buttons

      ############################################################################ constructing hgnc_filt using informations in the marker
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
      ############################################################################ end of constructing hgnc_filt dataframe of genename,id

      ### KEGG panel
      output$KEGG <- renderPrint({
        geneids<-hgnc_filt[rownames(toptable_sample[1:input$nogenes,]),2]
        keggtable<-kegga(geneids)
        kegg_ord<-keggtable[order(keggtable$P.DE),]
        kegg_ord[1:10,]
      })

      ### Markers panel
      output$Markers <- renderPrint({
        degenes<-rownames(toptable_sample[1:input$nogenes,])
        for (i in 1:nrow(hypergeom)) {
          hypergeom[i,1]<-names(markers_genelists_list)[i]
          hypergeom[i,2]<-phyper(length(intersect(degenes,markers_genelists_list[[i]])),length(markers_genelists_list[[i]]),ncol(scdata_subset)-1-length(markers_genelists_list[[i]]),length(degenes),lower.tail = F)
        }
        hypergeom_ord<-hypergeom[order(hypergeom$pvals),]
        hypergeom_ord[1:10,]
      })

      ### Msigdb panel
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



  #################################################################BOTTOM OF MAIN PANEL:
  #############################Adding tabset panel corresponds to each cluster (at the bottom of the main panel)
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
  ### end of tabset panel implementation

  ### update tabset panel according to the point selected
  observeEvent(event_data("plotly_click",priority = "event"), {
    cldat <- event_data("plotly_click")
    selectedi="No selection"
    selectedi <- filter(df,round(x1,5)==round(cldat$x,5))$y
    updateTabsetPanel(session, "switcher", selected = as.character(selectedi))
    i=df[,3][cldat$pointNumber]
    # if (previous_i!=-1){
    #   for (j in 1:length(c_intersections[[i]])){
    #     removeUI(
    #       selector = paste("#",c_intersections[[i]][j],"")
    #     )
    #   }
    # }
  })
  ############################# end of tabset panel

  ##########################################################################  Adding intersection buttons to corresponding tab panel
  #get the marker list intersections
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
          ) %>% layout(dragmode = "lasso",title=paste("Value of ",total_intersections[rr],sep=""))
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
  ###########################################################################  end of adding and maintaining intersection buttons
  ###################################################################################### END OF MAIN PANEL

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
