plots <- function(id, label='plots') {
ns = NS(id)

tabsetPanel(
    type = "tabs",
    id = ns("tabset"),
   
    
    tabPanel(
        "Main Plot",
        #uiOutput(ns("plots")),
        #splitLayout(
        #cellWidths = c('50%','50%'),
            jqui_resizable(plotlyOutput(ns('plot'),height = '800px',width = '1200px')),
            jqui_draggable(jqui_resizable(plotlyOutput(ns('plot2'),height='300px',width='500px')))    
        #)
        ,
        #tags$style(HTML(".tabbable > .nav > li > a  {background-color: aqua;  color:black}  ")),
        jqui_draggable(
            #tags$style(HTML(".tabbable > .nav > li > a  {background-color: aqua;  color:black}  ")),
        conditionalPanel(
            "output.plot",
            ns = ns,
            div(
                class = "cell_names_div",
                list(
                    actionButton(
                        ns("split_plot"),
                        "Split"
                    ),
                    actionButton(
                        ns("collapse_cell_names"),
                        "View additional info",
                        class = "collapsebtn"
                    ),
                    splitLayout(
                        htmlOutput(ns("clustering_info")),
                        htmlOutput(ns("cell_names_outp"))
            )))
        ),options = list(axis = 'y'))
     )

    
    #,
    # tabPanel(
    #     "Another Main Plot",
    # 
    #     jqui_resizable(plotlyOutput(ns("plotg"))),
    #     jqui_resizable(jqui_draggable(plotlyOutput(ns("plotg2"))))    
    #     
    #     ,
    #     conditionalPanel(
    #         "output.plotg",
    #         ns = ns,
    #         div(
    #             class = "cell_names_div",
    #             list(
    #                 actionButton(
    #                     ns("split_plot2"),
    #                     "Split"
    #                 )
    #                 ,
    #                 actionButton(
    #                     ns("collapse_cell_names2"),
    #                     "View additional info",
    #                     class = "collapsebtn"
    #                 )
    #                 ,
    #                 splitLayout(
    #                     htmlOutput(ns("clustering_info2")),
    #                     htmlOutput(ns("cell_names_outp2"))
    #                 )
    #             )
    #     )
    # ))

    
    )

}
