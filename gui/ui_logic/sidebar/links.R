github <- base64enc::dataURI(file="gui/ui_logic/icons/github.png", mime="image/png")
docs <- base64enc::dataURI(file="gui/ui_logic/icons/docs.png", mime="image/png")
mail <- base64enc::dataURI(file="gui/ui_logic/icons/mail.png", mime="image/png")

links = list(
    tags$li(a(href = 'https://github.com/ferrocactus/cellar',
            target = "_blank",
            icon("github"),
            title = "Github"),
            class = "dropdown"),
    tags$li(a(href = 'https://github.com/ferrocactus/cellar/blob/master/doc/cellar_guide.md',
            target = "_blank",
            icon("book"),
            title = "Documentation"),
            class = "dropdown"),
    tags$li(a(href = 'mailto:ehasanaj@cs.cmu.edu',
            target = "_blank",
            icon("envelope"),
            title = "Contact"),
            class = "dropdown")
)
