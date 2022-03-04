


#' @title View all plots for descriptive statistics.
#' 
#' @description 
#' 
#' This shiny module proposes an interface to choose and view the
#' module plots.
#' The list of these modules is available with the function [listPlotsModules()]
#' 
#' xxxxxx
#' 
#' @name all_plots
#' 
#' @examples 
#' if(interactive()){
#' library(QFeatures)
#' data(ft_na)
#' ui <- mod_all_plots_ui('allPlot')
#' 
#' server <- function(input, output, session) {
#'  data(ft_na)
#'  
#'  mod_all_plots_server('allPlot',
#'                       object = reactive({ft_na})
#'                       )
#'  }
#' shinyApp(ui=ui, server=server)
#' }
NULL




#' @rdname all_plots
#' @export
listPlotModules <- function(){
  ll <- ls('package:DaparToolshed')
  ll <- ll[grep('mod_plots_', ll)]
  ll <- gsub('_server', '', ll)
  ll <- gsub('_ui', '', ll)
  ll <- unique(ll)
  ll
}




#' @param id A `character(1)` for the 'id' of the shiny module. It must be 
#' the same as for the server function.
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs 
#' @rdname all_plots
#' @export
mod_all_plots_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    fluidPage(

      div( style="display:inline-block; vertical-align: middle; padding: 7px",
           uiOutput(ns('chooseDataset_ui'))
      ),
      div( style="display:inline-block; vertical-align: middle; padding: 7px",
           uiOutput(ns('ShowVignettes_ui'))
      ),

      br(),br(),br(),
      uiOutput(ns('ShowPlots_ui'))
  )
  )

}

#' @param id A `character(1)` for the 'id' of the shiny module. It must be 
#' the same as for the '*_ui' function.
#'
#' @param object A instance of the class `QFeatures`.
#' 
#' @importFrom base64enc dataURI
#' @importFrom shinyjs show hide hidden
#' 
#' @rdname all_plots
#' @export
#'
mod_all_plots_server <- function(id, object){

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    .width <- .height <- 40

    rv <- reactiveValues(
      current.plot = NULL,
      current.obj = NULL,
      current.indice = NULL,
      colData = NULL,
      conditions = NULL
    )

    
    observeEvent(input$rb, {
      tmp.list <- list.plots.module[-which(list.plots.module==input$rb)]
        shinyjs::show(paste0('div_', input$rb,'_large'))
        lapply(tmp.list, function(y){
          shinyjs::hide(paste0('div_', y,'_large'))
        })
      })



    output$ShowPlots_ui <- renderUI({
      lapply(list.plots.module, function(x){
        shinyjs::hidden(
          div(id=ns(paste0('div_', x, '_large')),
              do.call(paste0(x, '_ui'),
                      list(ns(paste0(x, '_large')))
                      )
              )
        )
      })
    })


    output$ShowVignettes_ui <- renderUI({

      tagList(
        tags$style(HTML("
            .shiny-input-radiogroup label {
                display: inline-block;
                text-align: center;
                margin-bottom: 20px;
            }
            .shiny-input-radiogroup label input[type='radio'] {
                display: inline-block;
                margin: 3em auto;
                margin-left: 10px;
            }
        ")),
        radioButtons(ns("rb"), "",
                     inline = TRUE,
                     choiceNames = lapply(list.plots.module, function(x){
                       img(src = dataURI(file = system.file('images', paste0(x, '.png'), package="MSPipelines"), mime="image/png"),
                           width='30px')
                     }),
                     choiceValues = list.plots.module,
                     selected = character(0)
                     )
      )
      })



    observeEvent(object(), {
      rv$colData <- SummarizedExperiment::colData(object())
      rv$metadata <- MultiAssayExperiment::metadata(object())
      rv$conditions <- SummarizedExperiment::colData(object())[['Condition']]
    })

    observeEvent(input$chooseDataset, {
      rv$current.indice <- which(names(object())==input$chooseDataset)
      rv$current.obj <- object()[[rv$current.indice]]
      })

    output$chooseDataset_UI <- renderUI({
      if (length(names(object())) == 0){
        choices <- list(' '=character(0))
      } else {
        choices <- names(object())
      }
      selectInput(ns('chooseDataset'), 'Dataset',
                  choices = choices,
                  selected = names(object())[length(object())],
                  width=200)
    })



    #Calls to server modules
    mod_plots_se_explorer_server('mod_plots_se_explorer_large',
                                 obj = reactive({rv$current.obj}),
                                 originOfValues = reactive({ rv$metadata[['OriginOfValues']] }),
                                 colData = reactive({ rv$colData })
    )


    mod_plots_intensity_server('mod_plots_intensity_large',
                               object=reactive({rv$current.obj}),
                               meta = reactive({ rv$metadata }),
                               conds = reactive({ rv$conditions }),
                               params = reactive({NULL}),
                               reset = reactive({FALSE}),
                               slave = reactive({FALSE}),
                               base_palette = reactive({
                                 DaparToolshed::Example_Palette(
                                   rv$conditions,
                                   DaparToolshed::Base_Palette(conditions = rv$conditions)
                                 )
                                 })
                               )


    mod_plots_pca_server('mod_plots_pca_large',
                         obj=reactive({rv$current.obj}),
                         coldata = reactive({ rv$colData })
    )


    mod_plots_var_dist_server('mod_plots_var_dist_large',
                              obj=reactive({rv$current.obj}),
                              conds = reactive({ rv$conditions}),
                              base_palette = reactive({
                                DaparToolshed::Example_Palette(
                                rv$conditions,
                                DaparToolshed::Base_Palette(conditions = rv$conditions)
                              )
                                })
                              )

    mod_plots_corr_matrix_server('mod_plots_corr_matrix_large',
                                 obj = reactive({rv$current.obj}),
                                 names = reactive({NULL}),
                                 gradientRate = reactive({0.9})
                                 )



    mod_plots_heatmap_server("mod_plots_heatmap_large",
                             obj = reactive({rv$current.obj}),
                             conds = reactive({ rv$conditions })
    )


    mod_plots_group_mv_server("mod_plots_group_mv_large",
                              obj = reactive({rv$current.obj}),
                              conds = reactive({ rv$colData }),
                              base_palette = reactive({
                                DaparToolshed::Example_Palette(
                                  rv$conditions,
                                  DaparToolshed::Base_Palette(conditions = rv$conditions)
                                )
                              })
    )

  })

}

## To be copied in the UI
# mod_all_plots_ui("all_plots_ui_1")

## To be copied in the server
# callModule(mod_all_plots_server, "all_plots_ui_1")

