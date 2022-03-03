# Module UI

#' @title   mod_var_dist_plot_ui and mod_var_dist_plot_server
#'
#' @description  A shiny Module.
#' 
#' @name cv_plots
#'
#' @examples
#'
#' #------------------------------
#' # plot a single plot
#' #------------------------------
#' 
#' CVDist(qData, conds)
#'
#' #------------------------------
#' # A shiny module
#' #------------------------------
#' 
#' if(interactive()){
#' library(QFeatures)
#' data(ft)
#' ui <- mod_cv_plot_ui('plot')
#' 
#' server <- function(input, output, session) {
#'  data(ft)
#'  conds <- design(ft)$Condition
#'  
#'  mod_cv_plot_server('plot',
#'                     obj = reactive({ft[[1]]}),
#'                     conds = reactive({conds}),
#'                     pal.name = reactive({'Dark2'})
#'                     )
#'  }
#' shinyApp(ui=ui, server=server)
#' }
NULL


#' @param id shiny id
#' @importFrom shiny NS tagList
#' @rdname cv_plots
#' @export
mod_cv_plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    helpText("Display the condition-wise distributions of the log-intensity CV (Coefficient of Variation)
               of the protein/peptides."),
    helpText("For better visualization, it is possible to zoom in by click-and-drag."),
    highchartOutput(ns("viewDistCV"),width = 600, height = 600)
  )
}



#' @param id shiny id
#' @param obj xxx
#' @param conds A `character()` representings the condition for 
#' each sample of the [QFeatures] object.
#' @param base_palette xxx
#' 
#' @importFrom shiny NS tagList
#' @rdname cv_plots
#' @export
mod_cv_plot_server <- function(id,
                               obj,
                               conds,
                               pal.name = NULL){


  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observe({
      req(obj())
      stopifnot(inherits(obj(), "SummarizedExperiment"))
    })



    viewDistCV <- reactive({
      req(obj())

      isolate({
        varDist <- CVDist(assay(obj()), conds(), pal.name)
      })
      varDist
    })


    output$viewDistCV <- renderHighchart({
      withProgress(message = 'Making plot', value = 100, {
        viewDistCV()
      })
    })

  })


}
