library(shiny)

selectionmenu <- menuItem(
    "Selection & Labeling",
    id = "selectionbtn",
    icon = icon("object-group"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_color",
                list(
                    selectInput(
                        "color",
                        "Select colour value:",
                        "cluster"
                    )
                ),
                HTML('<hr class="line">')
            ),

            div(
                class = "div_step div_subsets",
                list(
                    splitLayout(
                        cellWidths = c("50%", "50%"),
                        textInput(
                            "newsubset",
                            "New Subset",
                            placeholder = "Enter subset name"
                        ),
                        actionButton(
                            "store_lasso",
                            "Add Subset",
                            class="secondcol"
                        )
                    )
                ),
                HTML('<hr class="line">')
            ),

            div(
                class = "div_step div_labels",
                list(
                    splitLayout(
                        textInput(
                            "newlabelbox",
                            "New label",
                            placeholder = "Enter label to add"
                        ),
                        actionButton(
                            "labeladd",
                            "Add Label",
                            class="secondcol"
                        )
                    ),

                    splitLayout(
                        selectInput(
                            "newlabels",
                            "Select labels",
                            choices = 0
                        ),
                        actionButton(
                            "labelupd",
                            "Update Subset Labels",
                            class="secondcol"
                        )
                    )
                )
            )
        )
    )
)
