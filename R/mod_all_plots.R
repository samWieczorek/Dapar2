<<<<<<< HEAD
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
=======



#' @title View all plots for descriptive statistics.
#' 
#' @description 
#' 
#' This shiny module proposes an interface to choose and view the
#' module plots.
#' The list of these modules is available with the function [listPlotsModules()]
#' 
#' xxxxxx
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
#' 
#' @return A plot
#' 
<<<<<<< HEAD
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
=======
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
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
#' }
NULL




<<<<<<< HEAD



#' @rdname descriptive-statistics
#' @export
listPlotModules <- function(){
  ll <- ls('package:DaparToolshed')
  ll <- ll[grep('mod_ds_', ll)]
=======
#' @rdname all_plots
#' @export
listPlotModules <- function(){
  ll <- ls('package:DaparToolshed')
  ll <- ll[grep('mod_plots_', ll)]
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
  ll <- gsub('_server', '', ll)
  ll <- gsub('_ui', '', ll)
  ll <- unique(ll)
  ll
}




<<<<<<< HEAD
#' @param id A `character(1)` for the 'id' of the shiny module. It must be
#' the same as for the server function.
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs
#' @rdname descriptive-statistics
=======
#' @param id A `character(1)` for the 'id' of the shiny module. It must be 
#' the same as for the server function.
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs 
#' @rdname all_plots
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
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
<<<<<<< HEAD
    )
=======
  )
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
  )
  
}

<<<<<<< HEAD
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
  
=======
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

>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    .width <- .height <- 40
<<<<<<< HEAD
    
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
=======

    rv <- reactiveValues(
      current.plot = NULL,
      current.obj = NULL,
      current.indice = NULL,
      colData = NULL,
      conditions = NULL
    )

    list.plots.module <- listPlotModules()
    
    observeEvent(input$vignettes_rb, {
      tmp.list <- list.plots.module[-which(list.plots.module==input$vignettes_rb)]
        shinyjs::show(paste0('div_', input$vignettes_rb,'_large'))
        lapply(tmp.list, function(y){
          shinyjs::hide(paste0('div_', y,'_large'))
        })
      })



    output$ShowPlots_ui <- renderUI({
      lapply(list.plots.module, function(x){
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
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
      

<<<<<<< HEAD
    
    observeEvent(input$chooseDataset, {
      current.se(object()[[input$chooseDataset]])
    })

    output$chooseDataset_ui <- renderUI({
      req(object())
=======

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
        radioButtons(ns("vignettes_rb"), "",
                     inline = TRUE,
                     choiceNames = lapply(list.plots.module, function(x){
                       img(src = dataURI(file = system.file('images', paste0(x, '.png'), 
                                                            package="DaparToolshed"), 
                                         mime="image/png"),
                           width='30px')
                     }),
                     choiceValues = list.plots.module,
                     selected = character(0)
                     )
      )
      })



    observeEvent(object(), {
      rv$colData <- design(object())
      rv$metadata <- metadata(object())
      rv$conditions <- design(object())$Condition
    })

    observeEvent(input$chooseDataset, {
      #rv$current.indice <- which(names(object())==input$chooseDataset)
      rv$current.obj <- object()[[input$chooseDataset]]
      })

    output$chooseDataset_ui <- renderUI({
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
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
<<<<<<< HEAD
    
    
    
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
=======



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
                         obj = reactive({rv$current.obj}),
                         conds = reactive({ rv$conditions })
    )


    mod_plots_var_dist_server('mod_plots_var_dist_large',
                              obj = reactive({rv$current.obj}),
                              conds = reactive({ rv$conditions})
                              )

    mod_plots_corr_matrix_server('mod_plots_corr_matrix_large',
                                 qData = reactive({rv$current.obj})
                                 )



    mod_plots_heatmap_server("mod_plots_heatmap_large",
                             obj = reactive({rv$current.obj}),
                             conds = reactive({ rv$conditions })
    )


    mod_mv_plots_server("mod_plots_group_mv_large",
                              obj = reactive({rv$current.obj}),
                              conds = reactive({ rv$conditions })
>>>>>>> 815a43dd1eb217bab9139ef2898c98312a7bd60a
    )
    
  })
  
}
