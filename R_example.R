library(reticulate)
# Path to the folder of the python virtual environment
use_virtualenv("~/code/py3.7", required=TRUE)
##############################################################################
# Run this to import Pipeline into R, couldnt figure another way of doing it
source_python('__init__.py')
##############################################################################
# Read the data and the gene ids
ids <- read.csv('datasets/spleen/spleen.csv', nrows = 1, header = FALSE)
X <- read.csv('datasets/spleen/spleen.csv', skip = 1, header = FALSE)
##############################################################################
pipe <- Pipeline(X, config = 'configs/config.ini', col_ids = ids)
##############################################################################
# Run pipeline
pipe$run()
##############################################################################
# Labels are stored in
pipe$labels
##############################################################################
# Marker info is stored in pipe$markers
# This is a python dictionary which is converted to a named list in R
pipe$markers
# To access its elements, e.g., the lvl1 type of cluster "2", do
pipe$markers$"2"$lvl1_type
##############################################################################
# To run differential expression on a single subset of the data
# just store indices of the subset and pass to get_markers_subset
indices1 <- seq(1, 100)
dict <- pipe$get_markers_subset(indices1)
# You can also pass two subsets of indices and will return dict
indices2 <- seq(500, 900)
dict <- pipe$get_markers_subset(indices1, indices2)
##############################################################################
# To run constrained clustering, pass to update a list of new labels (all)
# let's generate some random labels (need at least 2 positive labels)
new_labels <- c(1)
new_labels <- rep(c(new_labels), dim(X)[1])
# let's change some labels to 0
for (x in seq(1, 1000)) {
    new_labels[x] = 0
}
# and let's change some to -1 (a labels of -1 will not be considered
# by the algorithm in clustering)
for (x in seq(1000, 2000)) {
    new_labels[x] = -1
}
# now update, use code=200 for seeded kmeans
pipe$update(new_labels, code=200)
# after update, can access new labels and new markers from pipe as above