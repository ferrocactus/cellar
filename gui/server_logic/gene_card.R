library(reticulate)

source_python("__init__.py")

getPage <- function(genename) {
    url <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                 genename, sep="")
    return(browseURL(url))
}

gene_card <- function(input, output, session) {
    observeEvent(input$search, {
        getPage(input$searchgene)
    })
}
