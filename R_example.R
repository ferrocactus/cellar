library(reticulate)
use_virtualenv("~/code/py3.7", required=TRUE)
source_python('example.py')

data <- load_data("spleen")
X <- data[1]
ids <- data[2]