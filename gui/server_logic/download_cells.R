source_python("__init__.py")
download_cells <- function(input, output, session,setNames, 
                           labelList, adata) 
{
  
  
  observeEvent(input$cell_subset_download,{ 
    subsets=""
    if (identical(input$cell_subset_download,NULL)){
      return()
    }
    else if (input$cell_subset_download==""){
      return()
    }
    
    labels=get_cluster_label_list(adata())
    clusters=as.character(input$cell_subset_download)
    #print(labels)
    #print(clusters)
    try(
      subsets<- validate_cluster_list(labels,clusters)
      
    )
    if (subsets!=""){
      ids=c()
      for (i in 1:length(subsets)){
        id=(subsets[[i]])
        
        obs=py_to_r(adata()$obs)
        ids=c(ids,obs$n_genes[which(obs$labels == id)])
      }
      #print(ids)
      
      allkeys=as.character(py_to_r(adata()$obs$index$to_numpy()))
      #print(allkeys)
      alldata=data.frame(allkeys,as.character(py_to_r(get_label_names(adata()))))
      colnames(alldata)<-c("Cell IDs","Updated Labels")
      ouput_labeldat<-alldata[ids,]
      output$download_cells <- downloadHandler(
        filename = "subsets.csv",
        content = function(file) {
          write.csv(ouput_labeldat, file, row.names = FALSE)
        }
      )
    }
    
    
    
  })
  #   if (pipe() == 0) return()
  #   if (!pipe()$has('labels')) return()
  #   if (!pipe()$has('x')) return()
  #   
  #    
  #   if (identical(input$cell_subset_download,NULL)){
  #       return()
  #   }
  #   else if (input$cell_subset_download==""){
  #       return()
  #   }
  #   else{
  #       ids=c()
  #       subsets=input$cell_subset_download
  #       if (substr(subsets,2,2)=="-"){
  #           a=strsplit(subsets,"-")[[1]][[1]]
  #           a=as.numeric(a)
  #           if (substr(subsets,3,3==""))
  #           {
  #               b=length(setNames())-1
  #           }
  #           else{
  #              b=strsplit(subsets,"-")[[1]][[2]]
  #              b=as.numeric(b)
  #           }
  #           for (i in a:b){
  #               ids=c(ids,setPts()[[i]])
  #           }
  #       }
  #       else{
  #           subsets=strsplit(subsets,",")[[1]]
  #           print(subsets)
  #           print("subsets")
  #           for (i in 1:length(subsets)){
  #               
  #               id=(subsets[[i]])
  #               #print(id)
  #               #print(setNames())
  #               #if (id < length((pipe()$get_cluster_names())[[2]])){
  #               ids=c(ids,setPts()[[which(setNames() == paste("Cluster_",id,sep=""))]])  
  #               #}
  #           }
  #       }
  #       #print(ids)
  #       #print(length(ids))
  #       
  # 
  #       #setname=pipe()$get_cluster_name(as.numeric(id))
  #       #id <- pipe()$get_cluster_id(input$cell_subset_download)
  #       #idx <- which(setNames() == input$cell_subset_download)
  #       #idx <- which(setNames() == input$cell_subset_download)
  #       #print("here")
  #       allkeys=c(1:length(pipe()$labels))
  #       alldata=data.frame(allkeys,as.character(pipe()$get_label_names()))#pipe()$labels)
  #       colnames(alldata)<-c("Keys","Updates Labels")
  #       ouput_labeldat<-alldata[ids,]
  #       
  #       #setname=gsub(':', '', setname)
  #       
  #   }
  #   
  #   output$download_cells <- downloadHandler(
  #     filename = "subsets.csv",
  #     content = function(file) {
  #       write.csv(ouput_labeldat, file, row.names = FALSE)
  #     }
  #   )
  #   
  #   
  # })
  
}