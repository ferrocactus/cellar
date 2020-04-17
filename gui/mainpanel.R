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
    actionButton("subset1", "store subset"),
    actionButton("subset2", "store subset"),
    actionButton("DEsubsets", "DE genes subsets"),
    uiOutput("genecard"),
    htmlOutput("inc"),
    splitLayout(
        textInput("searchgene", "Search Gene card",
                  placeholder = "Enter gene"),
        actionButton("search", "Search Card", class="sidebtn")
    ),

    fileInput(
        "file1", "Choose CSV File",
        multiple = FALSE,
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain", ".csv")
    )
)
