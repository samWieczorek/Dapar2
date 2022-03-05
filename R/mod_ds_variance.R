#' @param id shiny id
#' @importFrom shiny NS tagList
#' @rdname descriptive-statistics
#' @export
mod_ds_variance_ui <- function(id){
  ns <- NS(id)
  tagList(
    helpText("Display the condition-wise distributions of the log-intensity CV (Coefficient of Variation)
               of the protein/peptides."),
    helpText("For better visualization, it is possible to zoom in by click-and-drag."),
    highchartOutput(ns("viewDistCV"),width = 600, height = 600)
  )
}



#' @param id shiny id
#' @param object xxx
#' @param conds A `character()` representings the condition for 
#' each sample of the [QFeatures] object.
#' @param base_palette xxx
#' 
#' @importFrom shiny NS tagList
#' @rdname descriptive-statistics
#' @export
mod_ds_variance_server <- function(id,
                                   object,
                                   conds,
                                   pal.name = NULL){


  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observe({
      req(object())
      stopifnot(inherits(object(), "SummarizedExperiment"))
    })



    viewDistCV <- reactive({
      req(object())

      isolate({
        varDist <- CVDist(assay(object()), 
                          conds(), 
                          pal.name)
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
