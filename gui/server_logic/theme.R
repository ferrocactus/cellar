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
                selector = 'switcher'
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
            addCssClass(
                class = "odd_dark",
                selector = '.odd'
            )
            addCssClass(
                class = "datatables_dark",
                selector = '.datatables'
            )
            addCssClass(
                class = "dataTables_white",
                selector = '.dataTables_length'
            )
            addCssClass(
                class = "dataTables_white",
                selector = '.dataTables_filter'
            )
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
                selector = 'switcher'
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
            removeCssClass(
                class = "datatables_dark",
                selector = '.datatables'
            )
            removeCssClass(
                class = "dataTables_white",
                selector = '.dataTables_length'
            )
            removeCssClass(
                class = "dataTables_white",
                selector = '.dataTables_filter'
            )
        }
    })
}
