#' @title   Shiny module for density plots.
#'
#' @description  
#' 
#' A shiny Module.
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
#' data(ft)
#' ui <- fluidPage(
#' mod_plots_density_ui('plot')
#' )
#' 
#' server <- function(input, output, session) {
#'  data(ft)
#'  conds <- design(ft)$Condition
#'  legend <- design(ft)[["Sample.name"]]
#'   
#'  mod_plots_density_server('plot',
#'                            obj = reactive({ft}),
#'                            conds = reactive({conds}),
#'                            legend = reactive({legend}),
#'                            base_palette = reactive({Example_Palette()})
#'   )
#' }
#' 
#' 
#' shinyApp(ui=ui, server=server)
#' @export
#' @importFrom shiny NS tagList
#' @rdname descriptive-statistics-plots
mod_plots_density_ui <- function(id){
  ns <- NS(id)
  highchartOutput(ns("Densityplot"))
}


#' @export
#' @rdname descriptive-statistics-plots
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

