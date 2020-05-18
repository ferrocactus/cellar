datasetpanel <- div(
    id = "datasetpanel",
    class = "panelbody",

    selectInput(
        "dataset",
        "Choose a dataset:",
        choices = list.files("datasets")
    ),

    fileInput(
        "file1", "Choose CSV/h5ad File",
        multiple = FALSE,
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain", ".csv", ".h5ad")
    )
)
