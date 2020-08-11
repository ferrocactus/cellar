theme <- function(input, output, session, retheme) {
    observeEvent(input$theme_mode, {
        retheme(1)
    })

    observe({
        if (retheme() < 1) return()
        if (input$theme_mode == 'dark_mode') {
            addCssClass(
                class = "tab_dark_mode",
                selector = '.nav-tabs > li > a'
            )
            addCssClass(
                class = "nav_tabs_dark_mode",
                selector = '.nav-tabs'
            )
            addCssClass(
                class = "body_dark_mode",
                selector = '.content-wrapper'
            )
            addCssClass(
                class = "btn_dark",
                selector = '#ns-DEbuttons .btn.btn-default'
            )
            addCssClass(
                class = "cell_names_dark",
                selector = '#ns-collapse_cell_names'
            )
            runjs('
                $("<style>").text(".even { color:white !important; background-color: var(--bluegray) !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text(".odd { color:white !important; background-color: var(--darkbluegray) !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text("table.dataTable thead .sorting, .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button { color:white !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text("td > .btn-default { color:white !important; background-color: var(--darkbluegray) !important; } ").appendTo("head");
            ')
            runjs('
                $("<style>").text("td > .btn-default:hover { color:black !important; background-color: white !important; }").appendTo("head");
            ')
            runjs('
                $("<style>").text(".datatables button, .datatables input, .datatables select, .datatables textarea { background-color: var(--darkbluegray) !important; }").appendTo("head");
            ')
            runjs('
                $("<style>").text(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, .dataTables_wrapper .dataTables_paginate { color:white !important;}").appendTo("head");
            ')
        } else if (input$theme_mode == 'light_mode') {
            removeCssClass(
                class = "tab_dark_mode",
                selector = '.nav-tabs > li > a'
            )
            removeCssClass(
                class = "nav_tabs_dark_mode",
                selector = '.nav-tabs'
            )
            removeCssClass(
                class = "body_dark_mode",
                selector = '.content-wrapper'
            )
            removeCssClass(
                class = "btn_dark",
                selector = '#ns-DEbuttons .btn.btn-default'
            )
            removeCssClass(
                class = "cell_names_dark",
                selector = '#ns-collapse_cell_names'
            )
            runjs('
                $("<style>").text(".even { color:rgb(51, 51, 51) !important; background-color: white !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text(".odd { color:rgb(51, 51, 51) !important; background-color: rgb(249, 249, 249) !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text("table.dataTable thead .sorting, .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button { color:rgb(51, 51, 51) !important;}").appendTo("head");
            ')
            runjs('
                $("<style>").text("td > .btn-default { color:rgb(51, 51, 51) !important; background-color: #f4f4f4  !important;} ").appendTo("head");
            ')
            runjs('
                $("<style>").text("td > .btn-default:hover { color:rgb(51, 51, 51) !important; background-color: white !important; }").appendTo("head");
            ')
            runjs('
                $("<style>").text(".datatables button, .datatables input, .datatables select, .datatables textarea { background-color: white !important; } ").appendTo("head");
            ')
            runjs('
                $("<style>").text(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, .dataTables_wrapper .dataTables_paginate { color:rgb(51, 51, 51) !important;}").appendTo("head");
            ')
        }
        retheme(0)
    })
}
