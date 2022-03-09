
#' @param id xxx
#' @export
#' @importFrom shiny NS tagList
#' @importFrom highcharter highchartOutput
#' @rdname descriptive-statistics
mod_ds_mv_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    tagList(
      fluidRow(
        column(width = 4, 
               highchartOutput(ns("histo_MV")), 
               height="600px"),
        column(width = 4, 
               highchartOutput(ns("histo_MV_per_lines"))),
        column(width = 4, 
               highchartOutput(ns("histo_MV_per_lines_per_conditions")))
        )
      )
    )
}

#' @param id xxx
#' @param object An instance of the class [SummarizedExperiment]
#' @param conds A `Character()`
#' @param pal.name xxx
#'
#' @export
#' @importFrom highcharter renderHighchart
#' 
#' @rdname descriptive-statistics
#'
mod_ds_mv_server <- function(id,
                             object, 
                             conds, 
                             pal.name = reactive({NULL})){

  moduleServer(id, function(input, output, session){
    ns <- session$ns


    observe({
      req(object())
      stopifnot (inherits(object(), "SummarizedExperiment"))
    })



    output$histo_MV <- renderHighchart({
      req(object())
      conds()
      pal.name()

      withProgress(message = 'Making plot', value = 100, {
        tmp <- mvHisto(qData = assay(object()),
                       conds = conds(),
                       pal.name = pal.name())
      })
      tmp
    })



    output$histo_MV_per_lines <- renderHighchart({
      req(object())
      conds()

      isolate({
        withProgress(message = 'Making plot', value = 100, {
          tmp <- mvPerLinesHisto(qData = assay(object()))
        })
      })
      tmp
    })




    output$histo_MV_per_lines_per_conditions <- renderHighchart({
      req(object())
      conds()
      pal.name()

      withProgress(message = 'Making plot', value = 100, {
        tmp <- mvPerLinesHistoPerCondition(qData = assay(object()),
                                           conds = conds(),
                                           pal.name = pal.name())
      })
      tmp
    })


  })


}

## To be copied in the UI
# mod_plots_group_mv_ui("plots_group_mv_ui_1")

## To be copied in the server
# callModule(mod_plots_group_mv_server, "plots_group_mv_ui_1")

