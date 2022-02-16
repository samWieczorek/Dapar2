#' @title   mod_plots_density_ui and mod_plots_density_server
#'
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#' @param obj xxx
#' @param conds xxx
#' @param legend xxx
#' @param base_palette xxx
#'
#' @keywords internal
#' 
#' @return NA
#' 
#' @examples 
#' 
#' xxxx
#'
NULL

#' @importFrom shiny NS tagList
#'@rdname descriptives-statistics-plots
#'
mod_plots_density_ui <- function(id){
  ns <- NS(id)
  tagList(
    highchartOutput(ns("Densityplot"))
  )
}

#' @export
#'
#' @keywords internal
#'
#' @importFrom SummarizedExperiment assay
#' @rdname descriptives-statistics-plots
#'
mod_plots_density_server <- function(id, obj, conds, legend = NULL, base_palette = NULL){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observe({
      req(obj())
      if (class(obj()) != "SummarizedExperiment") { return(NULL) }
    })


    output$Densityplot <- renderHighchart({
      req(obj())
      req(conds())

      tmp <- NULL
      isolate({

        withProgress(message = 'Making plot', value = 100, {
          tmp <- DaparToolshed::densityPlotD_HC(qData = assay(obj()),
                                         conds = conds(),
                                         legend = legend(),
                                         palette = DaparToolshed::Base_Palette(conditions = conds()))
        })
      })
      tmp
    })


  })


}

## To be copied in the UI
# mod_plots_density_ui("plots_density_ui_1")

## To be copied in the server
# callModule(mod_plots_density_server, "plots_density_ui_1")

