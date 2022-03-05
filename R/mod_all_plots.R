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
#' #----------------------------------------
#' # Launch a single shiny module
#' #----------------------------------------
#' 
#' if(interactive()){
#'  data(ft)
#'  ui <- mod_ds_density_ui('plot')
#' 
#'  server <- function(input, output, session) {
#'   conds <- design(ft)$Condition
#'  
#'   mod_ds_density_server('plot',
#'                         qData = reactive({assay(ft, 1)}),
#'                         conds = reactive({conds})
#'                         )
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
    
    rv <- reactiveValues(
      current.plot = NULL,
      current.se = NULL,
      current.indice = NULL,
      colData = NULL,
      conditions = NULL
    )
    
    list.plots.module <- listPlotModules()
    
    observeEvent(input$vignettes_rb, {
      tmp.list <- list.plots.module[-which(list.plots.module == input$vignettes_rb)]
      shinyjs::show(paste0('div_', input$vignettes_rb,'_large'))
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
        radioButtons(ns("vignettes_rb"), "",
                     inline = TRUE,
                     choiceNames = lapply(list.plots.module, function(x){
                       fname <- gsub('mod_', '', x)
                       img(src = base64enc::dataURI(
                         file = system.file('images', 
                                            paste0(fname, '.png'),
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
      rv$current.se <- object()[[input$chooseDataset]]
    })
    
    output$chooseDataset_ui <- renderUI({
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
    
    
    
    #Calls to server modules
    mod_ds_explorer_server('mod_ds_explorer_large',
                           object = reactive({object()})
                           )
    
    
    mod_ds_intensity_server('mod_ds_intensity_large',
                            object = reactive({rv$current.se}),
                            exp.design = reactive({design(object()) })
                            )
    
    
    mod_ds_pca_server('mod_ds_pca_large',
                      object = reactive({rv$current.se}),
                      conds = reactive({design(object())$Condition})
                      )
    
    
    mod_ds_variance_server('mod_ds_variance_large',
                            object = reactive({rv$current.se}),
                            conds = reactive({design(object())$Condition})
                           )
    
    mod_ds_corrmatrix_server(id = 'mod_ds_corrmatrix_large',
                             object = reactive({rv$current.se})
                             )
    
    
    
    mod_ds_heatmap_server("mod_plots_heatmap_large",
                             object = reactive({rv$current.se}),
                             conds = reactive({design(object())$Condition})
    )
    
    
    mod_ds_mv_server("mod_ds_mv_large",
                        object = reactive({rv$current.se}),
                        conds = reactive({design(object())$Condition})
    )
    
  })
  
}

## To be copied in the UI
# mod_all_plots_ui("all_plots_ui_1")

## To be copied in the server
# callModule(mod_all_plots_server, "all_plots_ui_1")

