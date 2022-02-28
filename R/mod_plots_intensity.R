#' @title   mod_plots_intensity_plots_ui and mod_plots_intensity_plots_server
#'
#' @description  
#' 
#' A shiny Module.
#'
#' @name intensity_plots
#'
#' @examples 
#' 
#' 
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
               highchartOutput(ns("BoxPlot")),
               hidden(imageOutput(ns("violin_ui")))
      ),
      tags$div(style="display:inline-block; vertical-align: middle;",
               selectInput(ns("choosePlot"), "Choose plot",
                           choices=setNames(nm = c("violinplot", "boxplot")),
                           width='100px'),
               uiOutput(ns('slave_tracking_ui'))
      )
    )
  )
}

#' @param id xxx
#' @param object xxx
#' @param meta xxx
#' @param conds A `character()` of the names of conditions 
#' for each sample of the dataset.
#' @param base_palette xxx
#' @param params xxx
#' @param reset xxx
#' @param slave Default is FALSE.
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
                                       meta,
                                       conds,
                                       pal.name = NULL,
                                       params = NULL,
                                       reset = NULL,
                                       slave = FALSE){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    rv <- reactiveValues(
      indices = NULL,
      varTrack = NULL
    )



    rv$varTrack <- mod_tracking_server("slave_tracking",
                                       obj = reactive({object()}),
                                       keyId = reactive({meta()[['keyId']]}),
                                       params = reactive({params()}),
                                       reset = reactive({reset()}),
                                       slave = reactive({slave()})
    )

    output$slave_tracking_ui <- renderUI({
      req(!slave())
      req(typeDataset(object()) == 'protein')
      
      mod_tracking_ui(ns('slave_tracking'))
    })



    observeEvent(c(slave(), rv$varTrack()), ignoreInit = TRUE, ignoreNULL = FALSE, {
      if (slave()){
        switch(params()$typeSelect,
               ProteinList = rv$indices <- params()$list.indices,
               Random = rv$indices <- params()$rand.indices,
               Column = rv$indices <- params()$col.indices,
               None = rv$indices <- NULL
        )
      } else {
        tmp <- if (is.null(rv$varTrack()$typeSelect)) 
          'None'
        else 
          rv$varTrack()$typeSelect
        
        switch(tmp,
               ProteinList = rv$indices <- rv$varTrack()$list.indices,
               Random = rv$indices <- rv$varTrack()$rand.indices,
               Column = rv$indices <- rv$varTrack()$col.indices,
               None = rv$indices <- NULL
        )
        }
    })





    observeEvent(input$choosePlot, {
      shinyjs::toggle('viewViolinPlot', 
                      condition = input$choosePlot == 'violinplot')
      shinyjs::toggle('BoxPlot', 
                      condition = input$choosePlot == 'boxplot')
    })



    output$BoxPlot <- renderHighchart({
      object()
      rv$indices
      tmp <- NULL

      pattern <- paste0('test',".boxplot")
      withProgress(message = 'Making plot', value = 100, {
        tmp <- boxPlotD(object = object(),
                        exp.design = exp.design(),
                        keyId = rowData(object())[[ meta()[['keyId']] ]],
                        palette = DaparToolshed::Base_Palette(conditions=conds()),
                        subset.view = rv$indices)
      })
      tmp
    })


    output$viewViolinPlot<- renderImage({
      object()
      rv$indices
      tmp <- NULL

      # A temp file to save the output. It will be deleted after renderImage
      # sends it, because deleteFile=TRUE.
      outfile <- tempfile(fileext='.png')
      # Generate a png
      withProgress(message = 'Making plot', value = 100, {
        png(outfile)
        pattern <- paste0('test',".violinplot")
        tmp <- violinPlotD(assay(object()),
                           keyId = rowData(object())[[ meta()[['keyId']] ]],
                           conds = conds(),
                           palette = DaparToolshed::Base_Palette(conditions=conds()),
                           subset.view =  rv$indices)
        #future(createPNGFromWidget(tmp,pattern))
        dev.off()
      })
      tmp

      # Return a list
      list(src = outfile,
           alt = "This is alternate text")
    }, deleteFile = TRUE)


    return(reactive({rv$varTrack()}))


  })


}

## To be copied in the UI
# mod_plots_intensity_plots_ui("plots_intensity_plots_ui_1")

## To be copied in the server
# callModule(mod_plots_intensity_plots_server, "plots_intensity_plots_ui_1")

