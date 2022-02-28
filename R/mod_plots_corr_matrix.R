#' @title   mod_plots_group_mv_ui and mod_plots_group_mv_server
#'
#' @description  A shiny Module.
#' 
#' @name corrmatrix_plot
#'
#' @examples 
#' if (interactive()){
#' library(QFeatures)
#' library(DaparToolshed)
#' data(ft)
#' ui <- mod_corrmatrix_plot_ui('plot')
#' 
#' server <- function(input, output, session) {
#'  mod_corrmatrix_plot_server('plot',
#'                             qData = reactive({assay(ft[[1]])})
#'                             )
#'  }
#' shinyApp(ui=ui, server=server)
#' }
NULL


#' @param id xxx
#' @export
#' @importFrom shiny NS tagList
#' @importFrom shinyWidgets dropdownButton
#' @rdname corrmatrix_plot
mod_corrmatrix_plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('showValues_ui')),
    uiOutput(ns('rate_ui')),
    highchartOutput(ns("corrMatrix"), width = '600px',height = '500px')
  )
}

#' @param id xxx
#' @param qData xxx
#' @param rate xxx. Default value is 0.9
#' @param showValues Default is FALSE.
#'
#' @export
#' @rdname corrmatrix_plot
#'
mod_corrmatrix_plot_server <- function(id,
                                       qData, 
                                       rate = reactive({0.5}),
                                       showValues = reactive({FALSE})){

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      req(qData())
      stopifnot (inherits(qData(), "matrix"))
    })

    rv.corr <- reactiveValues(
      rate = NULL,
      showValues = FALSE
    )


    output$rate_ui <- renderUI({
      req(rv.corr$rate)
      sliderInput(ns("rate"),
                  "Tune to modify the color gradient",
                  min = 0,
                  max = 1,
                  value = rv.corr$rate,
                  step=0.01)
    })


    output$showValues_ui <- renderUI({
      checkboxInput(ns('showLabels'), 
                    'Show labels', 
                    value = rv.corr$showValues)

    })


    observeEvent(req(!is.null(input$showLabels)),{
      rv.corr$showValues <- input$showLabels
    })

    observeEvent(req(rate()),{
      rv.corr$rate <- rate()
    })

    observeEvent(req(input$rate),{
      rv.corr$rate <- input$rate
    })

    output$corrMatrix <- renderHighchart({
      req(qData())
      
      withProgress(message = 'Making plot', value = 100, {
        tmp <- corrMatrixPlot(qData = qData(),
                              rate = rv.corr$rate,
                              showValues = rv.corr$showValues)
      })
      
      tmp
    })


  })

}

## To be copied in the UI
# mod_plots_corr_matrix_ui("plots_corr_matrix_ui_1")

## To be copied in the server
# callModule(mod_plots_corr_matrix_server, "plots_corr_matrix_ui_1")
