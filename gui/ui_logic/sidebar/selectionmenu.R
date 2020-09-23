source("gui/ui_logic/sidebar/options.R")

selectionmenu <- function(id, label="selectionmenu") {
ns = NS(id)
menuItem(
    "Tools",
    id = "selectionbtn",
    icon = icon("object-group"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_subsets",
                list(
                    splitLayout(
                        #cellWidths = c("50%", "50%"),

                        textInput(
                            ns("newsubset"),
                            "New Subset",
                            placeholder = "Enter subset name"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "Store the selected points into a subset with this name.",
                                        placement= "bottom"
                                    )
                            ),

                        actionButton(
                            ns("store_lasso"),
                            "Add Subset",
                            class="secondcol"
                        )
                    )
                )
            ),

            div(
                class = "div_step div_labels",
                list(
                    splitLayout(
                        textInput(
                            ns("newlabelbox"),
                            "New label",
                            placeholder = "Enter label to add"
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip(
                                        "Enter the new label name you want to add",
                                        placement= "bottom"
                                    )
                            ),

                        actionButton(
                            ns("labeladd"),
                            "Add Label",
                            class = "secondcol"
                        )
                    ),

                    selectInput(
                        ns("tissue"),
                        "Select tissue",
                        choices = options$tissues
                    ) %>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Select the tissue of the cell you want to label",
                                    placement= "bottom"
                                )
                        ),

                    selectInput(
                        ns("newlabels"),
                        "Select cell type",
                        choices = c('None')
                    ) %>%
                        shinyInput_label_embed(
                            shiny::icon("info-circle") %>%
                                bs_embed_tooltip(
                                    "Select the cell type you want to label",
                                    placement= "bottom"
                                )
                        ),

                    splitLayout(
                        #cellWidths = c("49%", "1%","50%"),


                        #uiOutput("subset1_upd")

                    #uiOutput("subset1_upd"),
                    selectInput(
                        ns("subset1_upd"),
                        "Choose Subset",
                        choices = c("None")
                    ),
                        # redundant tooltip
                        # %>%
                        #     shinyInput_label_embed(
                        #         shiny::icon("info-circle") %>%
                        #             bs_embed_tooltip(
                        #                 "Select the tissue of the cell you want to label",
                        #                 placement= "bottom"
                        #             )
                        #     ),



                        actionButton(
                            ns("labelupd"),
                            "Update Subset Labels",
                            class = "secondcol"
                        )
                    )
                )
            ),

            div(
                class = "div_step div_upload_labels",
                list(
                    splitLayout(
                        fileInput(
                            ns("upload_labels"),
                            "Upload Labels",
                            multiple = FALSE,
                            accept = c(".csv")
                        ) %>%
                            shinyInput_label_embed(
                                shiny::icon("info-circle") %>%
                                    bs_embed_tooltip( title = paste0("Upload a csv containing labels"), placement = "bottom")
                            ),

                        actionButton(
                            ns("store_uploaded_labels"),
                            "Store Labels",
                            class = "secondcol"
                        )
                    )
                )
            )
        )
    )
)}
