library(shiny)

mainmenu <- menuItem(
    "Menu",
    id = "menubtn",
    icon = icon("th"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            h2("Colors"),
            selectInput("color", "Select colour value:", "cluster"),
            
            sliderInput("nogenes", "Select number of genes",
                        min = 1, max = 500, value = 10),
            
            
            h2("Selection Methods"),
            
            h4("Cluster Selection"),
            splitLayout(
                #cellWidths = c("20%", "40%", "40%"),
                selectInput("cluforanalysis", "Choose cluster", choices = c(""))
                
            ),
            
            splitLayout(
                cellWidths = c("40%", "40%", "20%"),
                actionButton("cluster_store","Store to Subset 1",
                             class="highbtn scdbtn"),
                actionButton("cluster_store2","Store to Subset 2",
                             class="highbtn scdbtn")
            ),
            
            
            h4("Lasso Selection"), 
            
            splitLayout(
                cellWidths = c("40%", "40%", "20%"),
                actionButton("subset1", "Store to Subset 1", class = "highbtn scdbtn"),
                actionButton("subset2", "Store to Subset 2",
                             class = "highbtn scdbtn")
                #actionButton("DEsubsets", "DE", class = "sidebtn scdbtn")
            ),
            
            
            
            
            
            h2("Assigning Labels"),
            splitLayout(
                textInput("newlabelbox", "New label",
                          placeholder = "Enter label to add"),
                actionButton("labeladd", "Add Label", class="sidebtn highbtn scdbtn")
            ),
            
            splitLayout(
                selectInput("newlabels", "Select labels", choices = 0),
                actionButton("labelupd", "Update Subset 1 Labels",
                             class="sidebtn highbtn scdbtn")
            ),
            
            h2("DE Genes"),
            
            actionButton("getdegenes", "Get DE genes", class="sidebtn longbtn"),
            
            h2("Gene Card"),
            uiOutput("genecard"),
            htmlOutput("inc"),
            splitLayout(
                textInput("searchgene", "Search Gene card",
                          placeholder = "Enter gene"),
                actionButton("search", "Search Card",
                             class="sidebtn highbtn scdbtn")
            )
            
            # splitLayout(
            #     cellWidths = c("40%", "40%", "20%"),
            #     actionButton("subset1", "Store Subset 1", class = "sidebtn"),
            #     actionButton("subset2", "Store Subset 2",
            #                  class = "sidebtn scdbtn"),
            #     actionButton("DEsubsets", "DE", class = "sidebtn scdbtn")
            # ),
            
            # splitLayout(
            #     cellWidths = c("50%", "50%"),
            #     selectInput("chgcluster","Cluster to Rename",choice = 0),
            #     textInput("newcluster", "New Name", value = "",
            #               width = NULL, placeholder = NULL)
            # ),
            # actionButton("chg","Apply Change"),
            #actionButton("disable","Disable buttons"),
            
            
            
            # actionButton("runanalysis", "Run analysis for clusters",
            #              class="sidebtn longbtn")
        )
    )
)
