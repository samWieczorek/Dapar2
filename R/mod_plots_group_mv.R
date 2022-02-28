
#' @param id xxx
#' @export
#' @importFrom shiny NS tagList
#' @importFrom highcharter highchartOutput
#' @rdname mv_plots
mod_mv_plots_ui <- function(id){
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
#' @param obj An instance of the class [SummarizedExperiment]
#' @param conds A `Character()`
#' @param pal.name xxx
#'
#' @export
#' @importFrom highcharter renderHighchart
#' 
#' @rdname mv_plots
#'
mod_mv_plots_server <- function(id, obj, conds, pal.name){

  moduleServer(id, function(input, output, session){
    ns <- session$ns


    observe({
      req(obj())
      stopifnot (inherits(obj(), "SummarizedExperiment"))
    })



    output$histo_MV <- renderHighchart({
      req(obj())
      conds()
      pal.name()

      withProgress(message = 'Making plot', value = 100, {
        tmp <- mvHisto(assay(obj()),
                       conds = conds(),
                       pal.name = pal.name())
      })
      tmp
    })



    output$histo_MV_per_lines <- renderHighchart({
      req(obj())
      conds()

      isolate({
        withProgress(message = 'Making plot', value = 100, {
          tmp <- mvPerLinesHisto(qData = assay(obj()))
        })
      })
      tmp
    })




    output$histo_MV_per_lines_per_conditions <- renderHighchart({
      req(obj())
      conds()
      pal.name()

      withProgress(message = 'Making plot', value = 100, {
        tmp <- mvPerLinesHistoPerCondition(qData = assay(obj()),
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

