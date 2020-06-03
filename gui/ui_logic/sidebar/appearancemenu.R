appearancemenu <- function(id, label='appearancemenu') {
  ns = NS(id)
  menuItem(
    "Appearance",
    id = "appearancebtn",
    icon = icon("palette"),
    startExpanded = FALSE,
    menuSubItem(
      icon = NULL,
      list(
        radioButtons(
            ns("theme_mode"),
            "Select theme:",
            c(
                "Light Mode" = "light_mode",
                "Dark Mode" = "dark_mode"
            ),
            inline = TRUE
        ),

        sliderInput(
          ns("dot_size"),
          "Select dot size",
          min = 1, max = 30, value = 5
        ),

        radioButtons(
            ns("show_names"),
            "Show names in plot:",
            c(
                "Don't show names" = "dont_show_names",
                "Show names" = "show_names"
            ),
            inline = TRUE
        )

      )
    )
  )}
