#' @param id A `character(1)` which is the id of the shiny module.
#' 
#' @export
#' @importFrom shiny NS tagList
#' @rdname descriptive-statistics
mod_ds_density_ui <- function(id){
  ns <- NS(id)
  tagList(
    highchartOutput(ns("plot_ui"))
  )
}


#' @param id A `character(1)` which is the id of the shiny module.
#' @param object xxx
#' @param conds xxx
#' @param legend xxx
#' @param base_palette xxx
#' @export
#' @rdname descriptive-statistics
#' 
mod_ds_density_server <- function(id,
                                    object, 
                                    conds,
                                    pal.name = reactive({NULL})
                                    ){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observe({
      req(object())
      stopifnot (inherits(object(), "SummarizedExperiment"))
    })


    output$plot_ui <- renderHighchart({
      req(object())
      req(conds())
      tmp <- NULL
      isolate({

        withProgress(message = 'Making plot', value = 100, {
          tmp <- densityPlot(object = object(),
                             conds = conds(),
                             pal.name = pal.name()
                             )
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

