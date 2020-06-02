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
        )

      )
    )
  )}
