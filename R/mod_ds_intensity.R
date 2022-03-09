#' @param id xxx
#' 
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs hidden
#' @importFrom stats setNames
#'
#' @return NA
#' @export
#' @rdname descriptive-statistics
mod_ds_intensity_ui <- function(id){
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
               mod_tracker_ui(ns('tracker'))
      )
    )
  )
}



#' @param id xxx
#' @param object xxx
#' @param exp.design xxx
#' @param ... xxx
#' 
#' @rdname descriptive-statistics
#'
#' @export
#' @importFrom grDevices png dev.off
#' @importFrom shinyjs toggle hidden
#'
#'@return NA
#'
mod_ds_intensity_server <- function(id,
                                    object,
                                    exp.design,
                                    ...){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    rv <- reactiveValues(
      indices = NULL
      )

    rv$indices <- mod_tracker_server("tracker",
                                       object = reactive({object()})
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

