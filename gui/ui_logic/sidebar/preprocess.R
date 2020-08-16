preprocessmenu <- function(id, label="preprocessmenu") {
ns = NS(id)
menuItem(
    "Preprocessing",
    id = "preprocessingbtn",
    icon = icon("filter"),
    startExpanded = FALSE,
    menuSubItem(
        icon = NULL,
        list(
            div(
                class = "div_step div_preprocess",
                list(
                    selectInput(
                        ns("preprocess_method"),
                        "Preprocessing (do not use with Server Datasets)",
                        c("Defaults", "Manual")
                    ),
                    conditionalPanel(
                        condition = "input.preprocess_method == 'Manual'",
                        ns = ns,
                        tags$p("Coming soon...")
                    ),
                    actionButton(
                        ns("runpreprocessbtn"),
                        "Run",
                        class = "sidebtn longbtn")))
        )
    )
)}