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

dataModal <- function(code = 0) {
    # Show dialog box when user inputs dataset name
    modalDialog(
        textInput("dataset_fname", "Enter dataset name."),
        if (code == 1) {
            div(tags$b(
                "Dataset exists.", style = "color: red;"
            ))
        } else if (code == 2) {
            div(tags$b(
                "Name should not be empty.", style = "color: red;"
            ))
        },
        footer = tagList(actionButton("okdataset", "OK"))
    )
}

datasetExists <- function(dataset, path) {
    # Check if dataset exists in path
    files <- list.files(path)

    if (length(files) > 0) {
        for (i in 1:length(files)) {
            if (tools::file_path_sans_ext(files[i]) == dataset) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

writeDataset <- function(path, name) {
    # File upload
    withProgress(
        message = 'Please wait...',
        value = 0, {
            incProgress(1/2, detail = paste("Processing file"))
            upload_file(name, path)
            incProgress(2/2, detail = paste("Saving file"))
        }
    )
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
