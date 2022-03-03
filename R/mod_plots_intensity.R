#' @title Intensity family plots
#' 
#' @description 
#' 
#' These functions are plotting different views of quantitative data.
#' Two plots are available:
#' 
#' - [boxplotD()] xxxx,
#' - [violinPlotD()] xxxx.
#' 
#' @author Samuel Wieczorek, Anais Courtier, Enora Fremy
#' 
#' @name intensity_plots
#' 
#' 
#' @examples
#' library(QFeatures)
#' data(ft)
#' 
#' #----------------------------------------------
#' # View violin plot with and without subset view
#' #----------------------------------------------
#' 
#' violinPlotD(ft[[1]], design(ft))
#' violinPlotD(ft[[1]], design(ft), subset = c(3,4,6), pal.name = 'Dark2')
#' 
#' #----------------------------------------
#' # xxxx
#' #----------------------------------------
#' 
#' boxPlotD(ft[[1]], design(ft))
#' boxPlotD(ft[[1]], design(ft), subset = c(2, 3), pal.name = 'Dark2')
#' 
#' # ----------------------------------------------------
#' # A shiny module which can track entities on the plots
#' # ----------------------------------------------------
#' 
#' if (interactive()){
#' library(QFeatures)
#' library(shiny)
#' library(shinyjs)
#' library(highcharter)
#' library(DaparToolshed)
#' 
#' data(ft)
#' 
#' 
#' ui <- tagList(
#'   mod_intensity_plots_ui('iplots'),
#'   uiOutput('show')
#' )
#' 
#' server <- function(input, output, session) {
#'   rv <- reactiveValues(
#'     tmp = NULL
#'   )
#'   
#'   mod_intensity_plots_server(id = 'iplots',
#'                              object = reactive({ft[[1]]}),
#'                              exp.design = reactive({design(ft)}),
#'                              pal.name = 'Dark2'
#'                              )
#'   
#' }
#' shinyApp(ui = ui, server = server)
#' }
#' 
NULL


#' @param id xxx
#' 
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs hidden
#'
#' @return NA
#' @export
#' @rdname intensity_plots
mod_intensity_plots_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    tags$div(
      tags$div(style="display:inline-block; vertical-align: middle;",
               uiOutput(ns("box_ui")),
               uiOutput(ns("violin_ui"))
      ),
      tags$div(style="display:inline-block; vertical-align: middle;",
               selectInput(ns("choosePlot"), 
                           "Choose plot",
                           choices = setNames(nm = c("violin", "box")),
                           width = '100px'),
               mod_tracking_ui(ns('tracker'))
      )
    )
  )
}



#' @param id xxx
#' @param object xxx
#' @param exp.design xxx
#' @param ... xxx
#' 
#' @rdname intensity_plots
#'
#' @export
#' @importFrom grDevices png
#' @importFrom shinyjs toggle hidden
#'
#'@return NA
#'
mod_intensity_plots_server <- function(id,
                                       object,
                                       exp.design,
                                       ...){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    rv <- reactiveValues(
      indices = NULL
      )

    rv$indices <- mod_tracking_server("tracker",
                                       obj = reactive({object()})
                                       )

    output$box_ui <- renderUI({
      if (input$choosePlot == 'box')
        highchartOutput(ns('box'))
      else
        hidden(highchartOutput(ns('box')))
    })
    
    output$box <- renderHighchart({
       withProgress(message = 'Making plot', value = 100, {
        tmp <- boxPlotD(object = object(),
                        exp.design = exp.design(),
                        subset = rv$indices(),
                        ...)
      })

    })


    output$violin_ui <- renderUI({
      if (input$choosePlot == 'violin')
        imageOutput(ns('violin'))
      else
        hidden(imageOutput(ns('violin')))
    })
    
    output$violin <- renderImage({
      
      # A temp file to save the output. It will be deleted after renderImage
      # sends it, because deleteFile=TRUE.
      outfile <- tempfile(fileext='.png')
      # Generate a png
      withProgress(message = 'Making plot', value = 100, {
        png(outfile)
        pattern <- paste0('test',".violinplot")
        tmp <- violinPlotD(object(),
                           exp.design = exp.design(),
                           subset = rv$indices(),
                           ...)
        #future(createPNGFromWidget(tmp,pattern))
        dev.off()
      })
      tmp
      
      # Return a list
      list(src = outfile,
           alt = "This is alternate text")
    }, deleteFile = TRUE)


  })


}

## To be copied in the UI
# mod_plots_intensity_plots_ui("plots_intensity_plots_ui_1")

## To be copied in the server
# callModule(mod_plots_intensity_plots_server, "plots_intensity_plots_ui_1")

