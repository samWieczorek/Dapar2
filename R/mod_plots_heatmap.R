#' @param id shiny id
#' @export
#' @importFrom shiny NS tagList
#' @return NA
#' @rdname heatmap_plot
mod_heatmap_plot_ui <- function(id){
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
#' @rdname heatmap_plot
#'
mod_heatmap_plot_server <- function(id, 
                                    obj, 
                                    conds, 
                                    width = 900){

  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    observe({
      req(obj())
      stopifnot (inherits(obj(), "SummarizedExperiment"))
    })


    limitHeatmap <- 20000
    height <- paste0(2*width/3,"px")
    width <- paste0(width,"px")

    output$DS_PlotHeatmap <- renderUI({
      req(obj())
      if (nrow(assay(obj())) > limitHeatmap){
        tags$p("The dataset is too big to compute the heatmap in a reasonable time.")
      }else {
        tagList(
          plotOutput(ns('heatmap_ui'), width = width, height = height)
        )
      }
    })



    output$heatmap_ui <- renderPlot({
      req(input$linkage)
      req(input$distance)
      
      isolate({
        withProgress(message = 'Making plot', value = 100, {
          heatmapD(qData = assay(obj()),
                   conds = conds(),
                   distance= input$distance,
                   cluster = input$linkage)
        })
      })
    })

  })

}

## To be copied in the UI
# mod_plots_heatmap_ui("plots_heatmap_ui_1")

## To be copied in the server
# callModule(mod_plots_heatmap_server, "plots_heatmap_ui_1")

