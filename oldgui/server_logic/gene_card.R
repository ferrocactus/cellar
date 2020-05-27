## getpage function for genecard
getPage <- function(genename) {
    url <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                 genename, sep="")
    return(browseURL(url))
}

getOutpNames <- function() {
    #TODO get all outp names
    outp_names=c("")
    for (i in 1:length(markers)){
        outp_names=c(outp_names,markers[[as.character(i-1)]][['outp_names']])
    }
}
