library(reticulate)

source_python("__init__.py")

gene_card <- function(input, output, session) {
    runjs("document.getElementById('ns-search').onclick = function() {
       window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + document.getElementById('ns-searchgene').value, '_blank'); };")
}
