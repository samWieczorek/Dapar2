#' @param id A `character(1)` which is the id of the shiny module.
#' 
#' @export
#' @importFrom shiny NS tagList
#' @rdname descriptive-statistics
mod_ds_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
    div(
      div(
        style="display:inline-block; vertical-align: middle; padding-right: 20px;",
        selectInput(ns("distance"),"Distance",
                    choices = list("Euclidean" ="euclidean",
                                   "Manhattan"="manhattan",
                                   "Maximum" = "maximum",
                                   "Canberra" = "canberra",
                                   "Binary" = "binary",
                                   "Minkowski" = "minkowski"),
                    selected = "euclidean",
                    width="150px")
      ),
      div(
        style="display:inline-block; vertical-align: middle; padding-right: 20px;",
        selectInput(ns("linkage"),"Linkage",
                    choices=list("Complete" = "complete",
                                 "Average"="average",
                                 "Ward.D"="ward.D",
                                 "Ward.D2"="ward.D2",
                                 "Single" = "single",
                                 "Centroid" = "centroid",
                                 "Mcquitty" = "mcquitty",
                                 "Median" = "median"),
                    selected='complete',
                    width="150px")
      ),
      tags$hr(),
      uiOutput(ns("DS_PlotHeatmap"))
    )
  )
}



#' @param id xxx
#' @param obj xxx
#' @param conds xxx
#' @param width xxx
#' @export
#' @rdname descriptive-statistics
#'
mod_ds_heatmap_server <- function(id, 
                                  object,
                                  conds,
                                  width = 900){

  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    observe({
      req(object())
      stopifnot (inherits(object(), "SummarizedExperiment"))
    })


    limitHeatmap <- 20000
    height <- paste0(2 * width / 3, "px")
    width <- paste0(width, "px")

    output$DS_PlotHeatmap <- renderUI({
      req(object())
      if (nrow(assay(object())) > limitHeatmap){
        tags$p("The dataset is too big to compute the heatmap in a reasonable time.")
      }else {
        tagList(
          plotOutput(ns('heatmap_ui'), 
                     width = width, 
                     height = height
                     )
        )
      }
    })



    output$heatmap_ui <- renderPlot({
      req(input$linkage)
      req(input$distance)
      
      isolate({
        withProgress(message = 'Making plot', value = 100, {
          heatmapD(qData = assay(object()),
                   conds = conds(),
                   distance= input$distance,
                   cluster = input$linkage)
        })
      })
    })

  })

}

