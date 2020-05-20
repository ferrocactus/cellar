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
            incProgress(1/4, detail = paste("Processing file"))
            l<-upload(name,path)
            filename<-l[[1]]
            df<-l[[2]]
            incProgress(2/4, detail = paste("Processing file"))
            mkdir(filename)
            incProgress(3/4, detail = paste("Saving file"))
            write_data(name,path,filename,df)   # use python to write h5ad file

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

    marker_n<-integer(length = length(marker_list_names))
    marker_nde<-integer(length = length(marker_list_names))
    marker_pvals<-double(length = length(marker_list_names))

    hypergeom<-data.frame(marker_list_names,marker_n)
    hypergeom<-data.frame(hypergeom,marker_nde)
    hypergeom<-data.frame(hypergeom,marker_pvals)
    colnames(hypergeom)<-c("NAME","N","Intersect Length","P Value")

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
