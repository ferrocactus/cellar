mainpanel <- div(
    id = "mainpanel",
    class = "panelbody",
    # Input: Selector for the gene to be plotted
    selectInput("color", "Select colour value:", "cluster"),

    sliderInput("nogenes", "Select number of genes",
                min = 1, max = 500, value = 10),

    splitLayout(
        textInput("newlabelbox", "New label",
                  placeholder = "Enter label to add"),
        actionButton("labeladd", "Add", class="sidebtn")
    ),

    splitLayout(
        selectInput("newlabels", "Select labels", choices = 0),
        actionButton("labelupd", "Update Labels", class="sidebtn")
    ),

    actionButton("getdegenes", "Get DE genes", class="sidebtn"),
    uiOutput("genecard"),
    htmlOutput("inc"),
    splitLayout(
        textInput("searchgene", "Search Gene card",
                  placeholder = "Enter gene"),
        actionButton("search", "Search Card", class="sidebtn")
    ),

    splitLayout(
        cellWidths = c("40%", "40%", "20%"),
        actionButton("subset1", "Store Subset 1", class = "sidebtn"),
        actionButton("subset2", "Store Subset 2", class = "sidebtn"),
        actionButton("DEsubsets", "DE", class = "sidebtn")
    ),
    
    # splitLayout(
    #     cellWidths = c("50%", "50%"),
    #     selectInput("chgcluster","Cluster to Rename",choice = 0),
    #     textInput("newcluster", "New Name", value = "", width = NULL, placeholder = NULL)
    # ),
    # actionButton("chg","Apply Change"),
    #actionButton("disable","Disable buttons"),
    
     
    selectInput("cluforanalysis", "Select clusters", choices = c("")),
    actionButton("runanalysis", "Run analysis for clusters", class="sidebtn"),
 
    fileInput(
        "file1", "Choose CSV/h5ad File",
        multiple = FALSE,
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain", ".csv", ".h5ad")
    )
)
