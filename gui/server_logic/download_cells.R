download_cells <- function(input, output, session, adata) {
    observeEvent(input$cell_subset_download,{
        output$download_cells <- downloadHandler(
        filename = "subsets.csv",
        content = function(file) {
            if (py_to_r(is_active(adata())) == FALSE) return()
            if (py_to_r(has_key(adata(), 'obs', 'labels')) == FALSE) return()

            subsets = ""
            if (identical(input$cell_subset_download, NULL)) {
                return()
            }

            else if (input$cell_subset_download == "") {
                return()
            }

            labels = (get_cluster_label_list(adata())) # all cluster labels
            clusters = as.character(input$cell_subset_download) # clusters selected

            try(
              subsets <- py_to_r(validate_cluster_list(labels, clusters)) ##check if the selected clusters are valid
            )

            # print(py_to_r(get_subsets(adata())))
            # print(py_to_r(get_cluster_label_list(adata())))
            # print(py_to_r(get_cluster_name_list(adata())))
            if (subsets != "") {
                outfile=data.frame("Cell_ID"="", "Label"="") ## initiate download file with an empty first row
                for (i in 1:length(subsets)) {        # for any input subset
                    id = as.character(subsets[[i]])               # get cluster id
                    allnames=py_to_r(get_subsets(adata()))
                    alllabels=py_to_r(get_cluster_label_list(adata()))
                    assigned_names=py_to_r(get_cluster_name_list(adata()))   # get updated cluster names
                    subsetname=allnames[which(allnames==paste0("Cluster_",as.character(id)))]  # eg. cluster_0
                    updated_name=assigned_names[which(alllabels==as.numeric(id))]
                    indices=py_to_r(adata()$uns['subsets'][subsetname])       # get indices of cells in that subset
                    updated_names=rep(as.character(updated_name),times=length(indices))
                    cell_ids=as.character(py_to_r(adata()$obs$index$to_numpy()))
                    indices=indices+1  # index in R starts from 1
                    #print(length(updated_names))
                    #print(length(indices))
                    #print(length(cell_ids[indices]))
                    #print(indices)
                    df = data.frame("Cell_ID"=cell_ids[indices], "Label"=updated_names)  # dataframe of this cluster
                    outfile=rbind(outfile,df)  # add this cluster to the output file

                }
                outfile=outfile[-1,]
                rownames(outfile)=1:length(outfile[,1])
                #allkeys = as.character(py_to_r(adata()$obs$index$to_numpy())) ## cell ids, i.e. row names
                #alllabels = py_to_r(get_label_names(adata()))
                #print(allkeys[1])
                #alldata = data.frame(allkeys, alllabels)
                #colnames(alldata) <- c("Cell IDs", "Updated Labels")
                #ouput_labeldat <- alldata[ids,]
                write.csv(outfile, file, row.names = FALSE)
            }
        }
      )
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