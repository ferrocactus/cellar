source("gui/ui_logic/header/history.R") # history
source("gui/ui_logic/header/title.R") # title
source("gui/ui_logic/header/info.R") # title

source("gui/ui_logic/sidebar/datasetmenu.R") #datasetmenu
source("gui/ui_logic/sidebar/appearancemenu.R") #appearancemenu
source("gui/ui_logic/sidebar/configmenu.R") #configmenu
source("gui/ui_logic/sidebar/alignmenu.R") #downloadmenu
source("gui/ui_logic/sidebar/selectionmenu.R") #mainmenu
source("gui/ui_logic/sidebar/analysismenu.R") #analysismenu
source("gui/ui_logic/sidebar/downloadmenu.R") #downloadmenu

source("gui/ui_logic/body/plots.R") #plots
source("gui/ui_logic/body/cell_names.R") #cell_names
source("gui/ui_logic/body/analysis.R") #analysis

header <- dashboardHeader(
    titleWidth = 400,
    title = title
    #info(id="ns")
)

sidebar <- dashboardSidebar(
    useShinyjs(),
    width = 400,
    sidebarMenu(
        datasetmenu(id='ns'),
        configmenu(id='ns'),
        alignmenu(id='ns'),
        selectionmenu(id='ns'),
        analysismenu(id='ns'),
        appearancemenu(id='ns'),
        downloadmenu(id='ns')
    )
)

body <- dashboardBody(
    useShinyjs(),
    tags$head(includeCSS("gui/ui_logic/styles/style.css")),
    tags$head(includeHTML(("gui/ui_logic/header/google-analytics.html"))),
    tags$script(HTML('
    $(document).ready(function() {
    $("header").find("nav").append(\'<span class="maintitle"> Cellar </span>\');
    })')),
    #tags$script(src = "gui/ui_logic/sidebar/anchor.js"),
    history(id='ns'),
    plots(id="ns"),
    cell_names(id="ns"),
    analysis(id="ns"),
    div(class = "footer",
        includeHTML("gui/ui_logic/body/copyright.html")
    )
)

ui <- dashboardPage(header, sidebar, body, skin = 'green', title = 'Cellar')
