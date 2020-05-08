library(reticulate)

source_python('__init__.py')
source_python('src/utils/read.py')

## getpage function for genecard
getPage <- function(genename) {
    url <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                 genename, sep="")
    return(browseURL(url))
}

intersect <- function(x, y) {
    y <- as.vector(y)
    unique(y[match(as.vector(x), y, 0L)])
}

# File upload
# File upload
writeDataset <- function(path, name) {
    withProgress(
        message = 'Please wait...',
        value = 0,
        {
            #incProgress(1/3, detail = paste("Catching file"))
            # 
            #f <- read.csv(path,
            #               #header = input$header,
            #               #sep = input$sep,
            #               #quote = input$quote)
            #)
            #f<-load_data(name)[1]
            #incProgress(1/4, detail = paste("Processing file"))
            
            fname <- strsplit(name, ".", fixed = TRUE)[[1]][1]
            type <- strsplit(name, ".", fixed = TRUE)[[1]][2]
            
            
            #dir.create(paste(getwd(), "/datasets/tmp", sep=""))
            #write.csv(f, paste(getwd(), "/datasets/tmp/tmp.csv", sep = ""))
            
            
            
            
            # incProgress(2/3, detail = paste("Saving file"))
            # dir.create(paste(getwd(), "/datasets/",as.character(fname), sep=""))
             if (type=="csv"){
                 f <- read.csv(path,
                                #header = input$header,
                                #sep = input$sep,
                                #quote = input$quote)
                 )
                 dir.create(paste(getwd(), "/datasets/",fname,"/", sep=""))
                 write.csv(f, paste(getwd(), "/datasets/",as.character(fname),"/",as.character(fname),".csv",sep=""))
             }
             else{
                 if (type=="h5ad"){
                     write_data(name,path)
                 }
             }
            
            
            
        }
    )
    #p
}

# Get Hypergeom
getHypergeom <- function(path) {
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

    #create hypergeon object
    marker_list_names <- names(markers_genelists_list)
    pvals <- double(length = length(markers_genelists_list))
    hypergeom <- data.frame(marker_list_names, pvals)

    return(hypergeom)
}

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
