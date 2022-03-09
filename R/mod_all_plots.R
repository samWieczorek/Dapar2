#' @title Bar plot of missing values per lines using highcharter.
#' 
#' @description 
#' 
#' This method plots a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins).
#' 
#' @param qData A `data.frame` or `matrix` that contains the quantitative data.
#' 
#' @param conds A `character()` of condition name for each sample. The 
#' length of 'conds' must be equal to the number of columns of 'qData'.
#' 
#' @param pal.name  A `character(1)` which is the name of the palette 
#' (from the package [RColorBrewer]) to use.
#' 
#' 
#' @details 
#' 
#' - distribution of the missing values per line,
#' 
#' - a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions,
#' 
#' - Histogram of missing values.
#' 
#' 
#' - Variance : Builds a densityplot of the CV of entities in numeric matrix.
#' The CV is calculated for each condition present in the dataset
#' (see the slot \code{'Condition'} in the \code{colData()} DataFrame)
#' 
#' - Heatmap:
#' 
#' 
#' The function [heatmapD()] 
#' 
#' 
#' The function [] is inspired from the function 'heatmap.2' 
#' that displays a numeric matrix. For more information, please refer to the help 
#' of the heatmap.2 function.
#' 
#' 
#' @section Missing values:
#' 
#' #' - distribution of the missing values per line,
#' 
#' - a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions,
#' 
#' - Histogram of missing values.
#' 
#' @name descriptive-statistics
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' library(DaparToolshed)
#' data(ft)
#' data(ft_na)
#' 
#' qData <- assay(ft, 1)
#' conds <- design(ft)$Condition
#' pal <- 'Dark2'
#' 
#' data(ft)
#' qData <- assay(ft, 1)
#' conds <- design(ft)$Condition
#' densityPlot(qData, conds)
#' 
#' legend <- design(ft)$Sample.name
#' densityPlot(qData, conds, pal.name = 'Dark2')
#' 
#' CVDist(qData, conds)
#' 
#' #----------------------------------------------
#' # Plots a heatmap for generic quantitative data
#' #----------------------------------------------
#' 
#' heatmapD(assay(ft, 1), conds)
#' 
#' #----------------------------------------------
#' # Plots a heatmap for missing values visualization
#' #----------------------------------------------
#' 
#' 
#' #' #----------------------------------------
#' # Plots a histogram of missing values
#' #----------------------------------------
#' 
#' mvPerLinesHisto(qData)
#' 
#' #----------------------------------------
#' # Plots a histogram of missing values
#' #----------------------------------------
#' 
#' mvHisto(qData, conds, pal.name = pal)
#' 
#' #----------------------------------------
#' # Plots a histogram of missing values
#' #----------------------------------------
#' 
#' 
#' 
#' library(QFeatures)
#' data(ft)
#' mv.heatmap(assay(ft, 1))
#' 
#' #----------------------------------------
#' # Launch a single shiny module
#' The list of all plots module is avaiable
#' # with the function [listPlotModules()].
#' # In the example, replace 'FOO' by the name of the module
#' # and add the necessaray parameters
#' #----------------------------------------
#' 
#' if(interactive()){
#'  data(ft)
#'  ui <- FOO_ui('plot')
#' 
#'  server <- function(input, output, session) {
#'   conds <- design(ft)$Condition
#'  
#'   FOO_server('plot',
#'              # Add appropriate parameters
#'              )
#'   }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
#' 
#' 
#' #----------------------------------------
#' # Launch a the shiny module for all plots
#' # Here, an example with the density plot
#' #----------------------------------------
#' 
#' if(interactive()){
#'  data(ft)
#'  ui <- mod_all_plots_ui('plot')
#' 
#'  server <- function(input, output, session) {
#'   mod_all_plots_server('plot',
#'                        object = reactive({ft})
#'                        )
#'   }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
NULL







#' @rdname descriptive-statistics
#' @export
listPlotModules <- function(){
  ll <- ls('package:DaparToolshed')
  ll <- ll[grep('mod_ds_', ll)]
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
#' @rdname descriptive-statistics
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
#' @rdname descriptive-statistics
#' @export
#'
mod_all_plots_server <- function(id, object){
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    .width <- .height <- 40
    
    ll.mods <- listPlotModules()
    
    current.se <- reactiveVal()
    btns.history <- reactiveVal(rep(0, length(ll.mods)))
    
    observeEvent(GetVignettesBtns(), {
      clicked <- which(btns.history() != GetVignettesBtns())
      shinyjs::show(paste0('div_', ll.mods[clicked],'_large'))
      lapply(ll.mods[-clicked], function(y){
         shinyjs::hide(paste0('div_', y,'_large'))
       })
      btns.history(GetVignettesBtns())
      
    })
  
    GetVignettesBtns <- reactive({
      unlist(lapply(ll.mods, 
                    function(x) input[[x]])
      )
    })
    
    output$ShowPlots_ui <- renderUI({
      lapply(ll.mods, function(x){
        shinyjs::hidden(
          div(id = ns(paste0('div_', x, '_large')),
              do.call(paste0(x, '_ui'),
                      list(ns(paste0(x, '_large')))
              )
          )
        )
      })
    })
    
    
    output$ShowVignettes_ui <- renderUI({
      lapply(ll.mods, function(x){
        actionButton(ns(x), 
                     label = tagList(
                       p(gsub('mod_ds_', '', x)),
                       tags$img(src = base64enc::dataURI(
                             file = system.file('images',
                                                paste0(gsub('mod_', '', x), '.png'),
                                                package='DaparToolshed'),
                             mime='image/png'),
                             height = "50px")
                       ),
                      style = 'padding: 0px;
                      border: none;
                      background-size: cover;
                           background-position: center;'
                         )
      }
                       )
        })
      

    
    observeEvent(input$chooseDataset, {
      current.se(object()[[input$chooseDataset]])
    })

    output$chooseDataset_ui <- renderUI({
      req(object())
      if (length(names(object())) == 0){
        choices <- list(' ' = character(0))
      } else {
        choices <- names(object())
      }
      selectInput(ns('chooseDataset'), 'Dataset',
                  choices = choices,
                  selected = names(object())[length(object())],
                  width = 200)
    })
    
    
    
    #
    # Calls to server modules
    #
    mod_ds_explorer_server('mod_ds_explorer_large',
                           object = reactive({object()}),
                           which.assay = reactive({input$chooseDataset})
                           )
    
    
    mod_ds_intensity_server('mod_ds_intensity_large',
                            object = reactive({current.se()}),
                            exp.design = reactive({design(object()) })
                            )
    
    
    mod_ds_pca_server('mod_ds_pca_large',
                      object = reactive({current.se()}),
                      conds = reactive({design(object())$Condition})
                      )
    
    
    mod_ds_variance_server('mod_ds_variance_large',
                            object = reactive({current.se()}),
                            conds = reactive({design(object())$Condition})
                           )
    
    mod_ds_corrmatrix_server(id = 'mod_ds_corrmatrix_large',
                             object = reactive({current.se()})
                             )
    
    
    
    mod_ds_heatmap_server("mod_plots_heatmap_large",
                             object = reactive({current.se()}),
                             conds = reactive({design(object())$Condition})
    )
    
    
    mod_ds_mv_server("mod_ds_mv_large",
                     object = reactive({current.se()}),
                     conds = reactive({design(object())$Condition})
    )
    
  })
  
}
