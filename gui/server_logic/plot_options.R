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
            console.log(gd);
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
"

js.reset_plot_height = "
    var plot = document.getElementById('ns-plot');
    var plot_height = document.getElementById('ns-plot_height').value;
    var upd = {height: parseInt(plot_height)};
    Plotly.relayout(plot, upd);
"

js.reset_theme = "
    var plot = document.getElementById('ns-plot');
    var theme_mode = $(\"input[name=ns-theme_mode]:checked\").val();
    console.log(theme_mode);
    var upd;
    if (theme_mode === 'dark_mode') {
        console.log(\"yes\");
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
    console.log(upd);
    Plotly.relayout(plot, upd);
"