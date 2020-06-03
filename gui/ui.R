source("gui/ui_logic/header/history.R")

source("gui/ui_logic/sidebar/datasetmenu.R") #datasetmenu
source("gui/ui_logic/sidebar/appearancemenu.R") #appearancemenu
source("gui/ui_logic/sidebar/configmenu.R") #configmenu
source("gui/ui_logic/sidebar/selectionmenu.R") #mainmenu
source("gui/ui_logic/sidebar/analysismenu.R") #analysismenu
source("gui/ui_logic/sidebar/downloadmenu.R") #downloadmenu

source("gui/ui_logic/body/plots.R") #plots
source("gui/ui_logic/body/cell_names.R") #cell_names
source("gui/ui_logic/body/analysis.R") #analysis

b64 <- base64enc::dataURI(file="gui/ui_logic/header/logo.png", mime="image/png")

header <- dashboardHeader(
    titleWidth = 400,
    title = div(list(
                img(src=b64, height='50px', style="float: left;"),
                p("Bar-Joseph Group", class = "title1"),
                tags$br(),
                p("School of Computer Science", class = "title2"),
                tags$br(),
                p("Carnegie Mellon University", class = "title3"),
                p("Cellar", class = "maintitle")
            ))
)

sidebar <- dashboardSidebar(
    useShinyjs(),
    width = 400,
    sidebarMenu(
        datasetmenu(id='ns'),
        configmenu(id='ns'),
        selectionmenu(id='ns'),
        analysismenu(id='ns'),
        appearancemenu(id='ns'),
        downloadmenu(id='ns')
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("gui/ui_logic/styles/style.css")),
    #tags$script(src = "gui/ui_logic/sidebar/anchor.js"),
    history(id='ns'),
    plots(id="ns"),
    cell_names(id="ns"),
    analysis(id="ns")
)

ui <- dashboardPage(header, sidebar, body, skin = 'green', title = 'Cellar')
