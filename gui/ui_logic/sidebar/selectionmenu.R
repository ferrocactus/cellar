library(shiny)

selectionmenu <- function(id, label="selectionmenu") {
ns = NS(id)
menuItem(
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
                        ns("color"),
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
                            ns("newsubset"),
                            "New Subset",
                            placeholder = "Enter subset name"
                        ),
                        actionButton(
                            ns("store_lasso"),
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
                            ns("newlabelbox"),
                            "New label",
                            placeholder = "Enter label to add"
                        ),
                        actionButton(
                            ns("labeladd"),
                            "Add Label",
                            class="secondcol"
                        )
                    ),

                    selectInput(
                        ns("tissue"),
                        "Select tissue",
                        choices = ''
                    ),

                    selectInput(
                        ns("newlabels"),
                        "Select cell type",
                        choices = ''
                    ),

                    splitLayout(
                        selectInput(
                            ns("subset1_upd"),
                            "Choose Subset",
                            choices = c("")
                        ),
                        actionButton(
                            ns("labelupd"),
                            "Update Subset Labels",
                            class="secondcol"
                        )
                    )
                )
            )
        )
    )
)}
