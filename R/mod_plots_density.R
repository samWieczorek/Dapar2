#' @title Shiny module for plots.
#'
#' @description  
#' 
#' A shiny Module.
#'
#' @name density-plots
#' 
#' @return A plot
#' 
#' @examples 
#' library(QFeatures)
#' data(ft)
#' ui <- mod_plots_density_ui('plot')
#' 
#' server <- function(input, output, session) {
#'  data(ft)
#'  conds <- design(ft)$Condition
#'  
#'  mod_plots_density_server('plot',
#'                           obj = reactive({ft}),
#'                           conds = reactive({conds}),
#'                           pal = reactive({'Dark1'})
#'                           )
#'  }
#' shinyApp(ui=ui, server=server)
NULL

#' @param id A `character(1)` which is the id of the shiny module.
#' 
#' @export
#' @importFrom shiny NS tagList
#' @rdname density-plots
mod_plots_density_ui <- function(id){
  fluidPage(
    ns <- NS(id),
    highchartOutput(ns("Densityplot"))
  )
}


#' @param id shiny id
#' @param obj xxx
#' @param conds xxx
#' @param legend xxx
#' @param base_palette xxx
#' @export
#' @rdname density-plots
#' 
mod_plots_density_server <- function(id, 
                                     obj, 
                                     conds, 
                                     legend = NULL, 
                                     base_palette = NULL){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observe({
      req(obj())
      stopifnot(inherits(obj(), "SummarizedExperiment"))
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

