library(reticulate)

source_python("__init__.py")

getPage <- function(genename) {
    url <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                 genename, sep="")
    return(browseURL(url, browser='firefox'))
}

gene_card <- function(input, output, session) {
    observeEvent(input$search, {
        c = Con()

        if (c$check_name(input$searchgene) == 1) {
            getPage(input$searchgene)
        } else {
            showNotification("Gene name does not exist.")
        }
    })
}
