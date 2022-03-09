#' @param id shiny id
#' @rdname descriptive-statistics
#' @export
#' @importFrom shiny NS tagList plotOutput
#' @importFrom highcharter highchartOutput
#' 
mod_plot_mvImputation_ui <- function(id){
  ns <- NS(id)
  tagList(
    tags$div(
      tags$div( style="display:inline-block; vertical-align: top; padding-right: 20px;",
                highchartOutput(ns("plot_viewNAbyMean"), width='600px')),
      tags$div( style="display:inline-block; vertical-align: top; padding-right: 20px;",
                plotOutput(ns("plot_showImageNA"), width='600px'))
    )
  )
}

#' @param id A `character(1)` which is the 'id' of the shiny module.
#' @param object An instance of a class `SummarizedExperiment`.
#' @param conds A `character()` of the name of conditions (one condition per sample).
#' @param title A `character(1)` which is
#' @param pal.name  A `character(1)` which is the name of the palette (from
#' the package [RColorBrewer] to use.
#' 
#' @rdname descriptive-statistics
#' @export
#' @importFrom highcharter renderHighchart
#' 
mod_plot_mvImputation_server <- function(id,
                                         object = reactive({NULL}),
                                         conds = reactive({NULL}),
                                         title = reactive({NULL}),
                                         pal.name = reactive({NULL})){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

  output$plot_viewNAbyMean <- renderHighchart({
      req(object())

      withProgress(message = 'Making MV Intensity plot', value = 100, {
        mvPlot.2(qData = assay(object()),
                 conds = conds(),
                 title = title(),
                 pal.name = pal.name())
      })
    })


    output$plot_showImageNA <- renderPlot({
      req(object())
      withProgress(message = 'Making MV Heatmap plot', value = 100, {
        mvImage(object = object(), 
                conds = conds())
      })

    })
  })
}

## To be copied in the UI
# mod_plots_mv_ui("plots_mv_ui_1")

## To be copied in the server
# callModule(mod_plots_mv_server, "plots_mv_ui_1")

