library(shiny)

mainmenu <- menuItem(
    "Menu",
    id = "menubtn",
    icon = icon("th"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            selectInput("color", "Select colour value:", "cluster"),

            sliderInput("nogenes", "Select number of genes",
                        min = 1, max = 500, value = 10),

            splitLayout(
                textInput("newlabelbox", "New label",
                          placeholder = "Enter label to add"),
                actionButton("labeladd", "Add", class="sidebtn highbtn scdbtn")
            ),

            splitLayout(
                selectInput("newlabels", "Select labels", choices = 0),
                actionButton("labelupd", "Update Labels",
                             class="sidebtn highbtn scdbtn")
            ),

            actionButton("getdegenes", "Get DE genes", class="sidebtn longbtn"),
            uiOutput("genecard"),
            htmlOutput("inc"),
            splitLayout(
                textInput("searchgene", "Search Gene card",
                          placeholder = "Enter gene"),
                actionButton("search", "Search Card",
                             class="sidebtn highbtn scdbtn")
            ),

            splitLayout(
                cellWidths = c("40%", "40%", "20%"),
                actionButton("subset1", "Store Subset 1", class = "sidebtn"),
                actionButton("subset2", "Store Subset 2",
                             class = "sidebtn scdbtn"),
                actionButton("DEsubsets", "DE", class = "sidebtn scdbtn")
            ),

            # splitLayout(
            #     cellWidths = c("50%", "50%"),
            #     selectInput("chgcluster","Cluster to Rename",choice = 0),
            #     textInput("newcluster", "New Name", value = "",
            #               width = NULL, placeholder = NULL)
            # ),
            # actionButton("chg","Apply Change"),
            #actionButton("disable","Disable buttons"),

            splitLayout(
                cellWidths = c("50%", "50%"),
                selectInput("cluforanalysis", "Choose cluster", choices = c("")),
                actionButton("cluster_label","Update cluster label",
                             class="highbtn scdbtn")
            ),

            actionButton("runanalysis", "Run analysis for clusters",
                         class="sidebtn longbtn")
        )
    )
)
