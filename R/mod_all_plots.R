

list.plots.module <- c(
  'mod_plots_se_explorer',
  'mod_plots_intensity',
  'mod_plots_pca',
  'mod_plots_var_dist',
  'mod_plots_corr_matrix',
  'mod_plots_heatmap',
  'mod_plots_group_mv'
  )


#' @title all_plots UI Function
#'
#' @description A shiny Module.
#' @export
#'
#' @param id xxx
#'
#' @importFrom shiny NS tagList
#'
mod_all_plots_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    fluidPage(

      div( style="display:inline-block; vertical-align: middle; padding: 7px",
           uiOutput(ns('chooseDataset_UI'))
      ),
      div( style="display:inline-block; vertical-align: middle; padding: 7px",
           uiOutput(ns('ShowVignettes'))
      ),

      br(),br(),br(),
      uiOutput(ns('ShowPlots'))
  )
  )

}

#' @title all_plots UI Function
#'
#' @description A shiny Module.
#' 
#' @export
#'
#' @param id xxxx
#' 
#' @param dataIn xxxx
#' 
#' @importFrom base64enc dataURI
#'
mod_all_plots_server <- function(id, dataIn){

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



    output$ShowPlots <- renderUI({
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


    output$ShowVignettes <- renderUI({

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
                     inline = T,
                     choiceNames = lapply(list.plots.module, function(x){
                       img(src = base64enc::dataURI(file=system.file('images', paste0(x, '.png'), package="MSPipelines"), mime="image/png"),
                           width='30px')
                     }),
                     choiceValues = list.plots.module,
                     selected = character(0)
                     )
      )
      })



    observeEvent(dataIn(), {
      rv$colData <- SummarizedExperiment::colData(dataIn())
      rv$metadata <- MultiAssayExperiment::metadata(dataIn())
      rv$conditions <- SummarizedExperiment::colData(dataIn())[['Condition']]
    })

    observeEvent(input$chooseDataset, {
      rv$current.indice <- which(names(dataIn())==input$chooseDataset)
      rv$current.obj <- dataIn()[[rv$current.indice]]
      })

    output$chooseDataset_UI <- renderUI({
      if (length(names(dataIn())) == 0){
        choices <- list(' '=character(0))
      } else {
        choices <- names(dataIn())
      }
      selectInput(ns('chooseDataset'), 'Dataset',
                  choices = choices,
                  selected = names(dataIn())[length(dataIn())],
                  width=200)
    })



    #Calls to server modules
    mod_plots_se_explorer_server('mod_plots_se_explorer_large',
                                 obj = reactive({rv$current.obj}),
                                 originOfValues = reactive({ rv$metadata[['OriginOfValues']] }),
                                 colData = reactive({ rv$colData })
    )


    mod_plots_intensity_server('mod_plots_intensity_large',
                               dataIn=reactive({rv$current.obj}),
                               meta = reactive({ rv$metadata }),
                               conds = reactive({ rv$conditions }),
                               params = reactive({NULL}),
                               reset = reactive({FALSE}),
                               slave = reactive({FALSE}),
                               base_palette = reactive({
                                 Dapar2::Example_Palette(
                                   rv$conditions,
                                   Dapar2::Base_Palette(conditions = rv$conditions)
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
                                Dapar2::Example_Palette(
                                rv$conditions,
                                Dapar2::Base_Palette(conditions = rv$conditions)
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
                                Dapar2::Example_Palette(
                                  rv$conditions,
                                  Dapar2::Base_Palette(conditions = rv$conditions)
                                )
                              })
    )

  })

}

## To be copied in the UI
# mod_all_plots_ui("all_plots_ui_1")

## To be copied in the server
# callModule(mod_all_plots_server, "all_plots_ui_1")

