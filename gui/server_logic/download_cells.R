download_cells <- function(input, output, session,setNames, setPts,
                          labelList, pipe) 
{
  
  
  
  observeEvent(input$cell_subset_download,{ 
    if (pipe() == 0) return()
    if (!pipe()$has('labels')) return()
    if (!pipe()$has('x')) return()
    
    
    if (input$cell_subset_download =="All"){
        ouput_labeldat<-data.frame(pipe()$x,pipe()$labels)
        colnames(ouput_labeldat)<-c(pipe()$col_ids,"Updates Labels")
        filename="All.csv"
        comp=NULL
        compname="None.csv"
    }
    else if (input$cell_subset_download == "None"){
        ouput_labeldat<-NULL
        filename="None.csv"
        compname="All.csv"
        comp<-data.frame(pipe()$x,pipe()$labels)
        colnames(comp)<-c(pipe()$col_ids,"Updates Labels")
    }
    else if (identical(input$cell_subset_download,NULL)){
      return()
    }
    else{
        id<-strsplit(input$cell_subset_download,"_")[[1]][[2]]
        id=as.numeric(id)
        setname=pipe()$get_cluster_name(as.numeric(id))
        #id <- pipe()$get_cluster_id(input$cell_subset_download)
        idx <- which(setNames() == input$cell_subset_download)

        #print(setNames())
        #print(idx)
        #print(setname)

        #idx <- which(setNames() == input$cell_subset_download)
        keys <- setPts()[[idx]]
        compkey=setdiff((1:length(pipe()$col_ids)),keys)
        alldata=data.frame(pipe()$x,pipe()$labels)
        colnames(alldata)<-c(pipe()$col_ids,"Updates Labels")
        ouput_labeldat<-alldata[keys,]
        
        setname=gsub(':', '', setname)
        filename=paste0(setname,".csv")
        compname=paste0(setname,"_complement.csv")
        comp=alldata[,compkey]
    }
    
    output$download_cells <- downloadHandler(
      filename = filename,
      content = function(file) {
        write.csv(ouput_labeldat, file, row.names = FALSE)
      }
    )
    
     output$download_comp <- downloadHandler(
       filename = compname,
       content = function(file) {
         write.csv(comp, file, row.names = FALSE)
       }
     )
    
    
  })
  
}