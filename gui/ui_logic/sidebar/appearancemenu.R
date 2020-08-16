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
        ) ,
        # %>%
        #   shinyInput_label_embed(
        #     shiny::icon("info-circle") %>%
        #       bs_embed_tooltip("Choose a theme of the app", placement="bottom",position="right",aligh="right")
        #   ),
        
        
        radioButtons(
          ns("show_names"),
          "Show names in plot:",
          c(
            "Don't show names" = "dont_show_names",
            "Show names" = "show_names"
          ),
          inline = TRUE
        ) ,
        # %>%
        #   shinyInput_label_embed(
        #     shiny::icon("info-circle") %>%
        #       bs_embed_tooltip(paste0("show/hide cluster names in the plot"), placement="bottom",position="right")
        #   ),
        
        sliderInput(
          ns("dot_size"),
          "Select dot size",
          min = 1, max = 30, value = 5
        ) ,
        # %>%
        #   shinyInput_label_embed(
        #     shiny::icon("info-circle") %>%
        #       bs_embed_tooltip(paste0("Select the size of the dots in the plot"), placement="bottom",position="right")
        #   ),
        
        sliderInput(
          ns("plot_height"),
          "Select plot height",
          min = 400, max = 800, value = 600
        ) 
        # %>%
        #   shinyInput_label_embed(
        #     shiny::icon("info-circle") %>%
        #       bs_embed_tooltip(paste0("Adjust the height of the main plot"), placement="bottom",position="right")
        #   )
        
      )
    )
  )}
