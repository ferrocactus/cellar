source("gui/ui_logic/header/history.R")
source("gui/ui_logic/header/title.R")
source("gui/ui_logic/header/info.R")

source("gui/ui_logic/sidebar/datasetmenu.R")
source("gui/ui_logic/sidebar/preprocess.R")
source("gui/ui_logic/sidebar/appearancemenu.R")
source("gui/ui_logic/sidebar/clusterimenu.R")
source("gui/ui_logic/sidebar/clusteriimenu.R")
source("gui/ui_logic/sidebar/alignmenu.R")
source("gui/ui_logic/sidebar/selectionmenu.R")
source("gui/ui_logic/sidebar/analysismenu.R")
source("gui/ui_logic/sidebar/downloadmenu.R")
source("gui/ui_logic/sidebar/links.R")

source("gui/ui_logic/body/plots.R")
source("gui/ui_logic/body/analysis.R")

source("gui/ui_logic/tooltips.R")

header <- dashboardHeader(
    titleWidth = 400,
    title = title,
    links[[1]],
    links[[2]],
    links[[3]]
)

sidebar <- dashboardSidebar(
    useShinyjs(),
    width = 400,
    sidebarMenu(
        datasetmenu(id='ns'),
        preprocessmenu(id='ns'),
        clusterimenu(id='ns'),
        clusteriimenu(id='ns'),
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
    tags$head(tags$link(rel = "icon", type = "image/x-icon",
              href = base64enc::dataURI(
                  file="gui/ui_logic/icons/favicon.ico", mime="image/ico"))),
    includeScript("gui/scripts/cellar.js"),
    includeScript("gui/scripts/anchor_scroll_fix.js"),
    tags$script(HTML('
        $("link[href*=\'_all-skins.min.css\']").remove();
    ')),
    tags$head(includeCSS("gui/ui_logic/styles/_all-skins.min.css")),
    history(id='ns'),
    plots(id="ns"),
    analysis(id="ns"),
    tooltips(id="ns"),
    div(class = "footer",
        includeHTML("gui/ui_logic/body/copyright.html")
    )
)

ui <- dashboardPage(header, sidebar, body, skin = 'purple', title = 'Cellar')
