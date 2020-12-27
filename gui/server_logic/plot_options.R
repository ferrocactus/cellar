all_symbols <- c(
    'star-triangle-down','circle','square',"diamond",
    "x",'triangle-up','triangle-down','hexagon','asterisk',
    'diamond-cross',"square-cross","circle-cross","circle-x",
    "star-square","star","star-triangle-up","star-square",
    "star-diamond","diamond-tall","diamond-wide","hourglass",
    'bowtie','pentagon','hexagram-dot','triangle-se','y-right',
    'hexagon2','octagon','triangle-nw','triangle-sw')

plot_options_btn <- list(
    name = "Options",
    #icon = icon("cog"),
    click = htmlwidgets::JS(
        "function(gd, upd) {
            var upd = {showlegend: false};
            Plotly.restyle(gd, upd);
        }"
    )
)

js.reset_marker_size = "
    var plot = document.getElementById('ns-plot');
    var marker_size = document.getElementById('ns-dot_size').value;
    var upd = {'marker.size': marker_size};
    Plotly.restyle(plot, upd);

    var plot2 = document.getElementById('ns-plot2');
    if (plot2 != null) {
        Plotly.restyle(plot2, upd);
    }

    var plots = document.getElementById('ns-tabset').children;
    for (var i = 1; i < plots.length; ++i) {
        var plot_num = plots[i].childNodes[1].innerText.split(' ')[1];
        plot = document.getElementById('ns-frozen_plot' + plot_num);
        Plotly.restyle(plot, upd);
    }
"

js.reset_plot_height = "
    var plot_height = document.getElementById('ns-plot_height').value;
    var plot = document.getElementById('ns-plot');
    var upd = {height: parseInt(plot_height)};
    Plotly.relayout(plot, upd);

    var plot2 = document.getElementById('ns-plot2');
    if (plot2 != null) {
        Plotly.relayout(plot2, upd);
    }

    var plots = document.getElementById('ns-tabset').children;
    for (var i = 1; i < plots.length; ++i) {
        var plot_num = plots[i].childNodes[1].innerText.split(' ')[1];
        plot = document.getElementById('ns-frozen_plot' + plot_num);
        Plotly.relayout(plot, upd);
    }
"

js.reset_theme = "
    var theme_mode = $(\"input[name=ns-theme_mode]:checked\").val();
    var upd;
    if (theme_mode === 'dark_mode') {
        upd = {
            plot_bgcolor: 'rgb(44, 59, 65)',
            paper_bgcolor: 'rgb(44, 59, 65)',
            'font.color': 'rgb(255, 255, 255)',
            'xaxis.gridcolor': 'rgb(85, 85, 85)',
            'xaxis.zerolinecolor': 'rgb(125, 125, 125)',
            'yaxis.gridcolor': 'rgb(85, 85, 85)',
            'yaxis.zerolinecolor': 'rgb(125, 125, 125)'
        };
    } else {
        upd = {
            plot_bgcolor: null,
            paper_bgcolor: null,
            'font.color': null,
            'xaxis.gridcolor': null,
            'xaxis.zerolinecolor': null,
            'yaxis.gridcolor': null,
            'yaxis.zerolinecolor': null
        };
    }
    var plot = document.getElementById('ns-plot');
    Plotly.relayout(plot, upd);

    var plot2 = document.getElementById('ns-plot2');
    if (plot2 != null) {
        Plotly.relayout(plot2, upd);
    }

    var plots = document.getElementById('ns-tabset').children;
    for (var i = 1; i < plots.length; ++i) {
        var plot_num = plots[i].childNodes[1].innerText.split(' ')[1];
        plot = document.getElementById('ns-frozen_plot' + plot_num);
        Plotly.relayout(plot, upd);
    }
"

theme_plot <- function(p, theme_mode) {
    if (theme_mode == 'dark_mode') {
        plot_bgcolor = 'rgb(44, 59, 65)'
        paper_bgcolor = 'rgb(44, 59, 65)'
        font = list(color = 'rgb(255, 255, 255)')
        xaxis = list(gridcolor = 'rgb(85, 85, 85)',
                     zerolinecolor = 'rgb(125, 125, 125)')
        yaxis = list(gridcolor = 'rgb(85, 85, 85)',
                     zerolinecolor = 'rgb(125, 125, 125)')
    } else {
        plot_bgcolor = 'rgb(255, 255, 255)'
        paper_bgcolor = 'rgb(255, 255, 255)'
        font = list(color = 'rgb(0, 0, 0)')
        xaxis = list(gridcolor = '#eee',
                     zerolinecolor = '#444')
        yaxis = list(gridcolor = '#eee',
                     zerolinecolor = '#444')
    }

    p <- p %>% layout(
        plot_bgcolor = plot_bgcolor,
        paper_bgcolor = paper_bgcolor,
        font = font,
        xaxis = xaxis,
        yaxis = yaxis)

    return(p)
}