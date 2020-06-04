b64 <- base64enc::dataURI(file="gui/ui_logic/header/logo.png", mime="image/png")

title = div(list(
    img(src=b64, height='50px', style="float: left;"),
    p("Bar-Joseph Group", class = "title1"),
    tags$br(),
    p("School of Computer Science", class = "title2"),
    tags$br(),
    p("Carnegie Mellon University", class = "title3")))
