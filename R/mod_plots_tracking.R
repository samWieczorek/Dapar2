#' @title  Tracking of entities within plots
#'
#' @description  This shiny module offers a UI to select a subset of a dataset
#' and superimpose quantitative values of this selection on the complete plot
#' Three modes of selection are implemented:
#'
#' - 'Protein list': xxx,
#' - 'Random': xxx,
#' - 'Specific column': xxx
#' 
#' @name tracking
#' @examples 
#' if (interactive()){
#' library(QFeatures)
#' library(DaparToolshed)
#' data(ft)
#' ui <- mod_tracking_ui('track')
#' 
#' server <- function(input, output, session) {
#'  mod_tracking_server(id = 'track',
#'                      object = reactive({ft[[1]]})
#'                      )
#'  }
#' shinyApp(ui, server)
#' }
NULL

 
#' @param id shiny id
#' @export
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs hidden
#' 
#' @rdname tracking
#' 
#' @return NA
#'
mod_tracking_ui <- function(id){
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    actionButton(ns('rst_btn'), 'Reset'),
    uiOutput(ns('typeSelect_ui')),
    hidden(uiOutput(ns('listeSelect_ui'))),
    hidden(uiOutput(ns('randSelect_ui'))),
    hidden(uiOutput(ns('colSelect_ui')))
  )
}

#' @param id xxx
#' @param object A instance of the class `SummarizedExperiment`
#' @param reset xxx
#' 
#' @rdname tracking
#'
#' @export
#'
#' @importFrom shinyjs toggle hidden show hide
#' 
#' @return A `list()` of integers
#'
mod_tracking_server <- function(id,
                                object){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    rv.track <- reactiveValues(
      typeSelect = "None",
      listSelect = NULL,
      randSelect = '',
      colSelect = NULL,
      indices = NULL)

    
    output$typeSelect_ui <- renderUI({
      tmp <- c("None", "ProteinList", "Random", "Column")
      nm <- c("None", "Protein list", "Random", "Specific Column")
      
      selectInput(ns("typeSelect"), 
                  "Type of selection",
                  choices = setNames(tmp, nm),
                  selected = rv.track$typeSelect,
                  width = '130px')
    })

    
    output$listSelect_ui <- renderUI({
      selectInput(ns("listSelect"),
                  "Select protein",
                  choices = c('None'),
                  multiple = TRUE,
                  selected = rv.track$listSelect,
                  )
    })
    
    output$colSelect_ui <- renderUI({
      selectInput(ns("colSelect"),
                  "Column of rowData",
                  choices = c(''),
                  selected = rv.track$colSelect
                  )
    })
    
    output$randSelect_ui <- renderUI({
      textInput(ns("randSelect"),  
                "Random", 
                value = rv.track$randSelect,
                width = ('120px'))
    })
    
    observeEvent(req(input$typeSelect),{
      rv.track$typeSelect <- input$typeSelect
      
      print(input$typeSelect)
      
      shinyjs::toggle(id = ns('listSelect'), 
                      condition = rv.track$typeSelect == 'ProteinList')
      shinyjs::toggle(id = 'randSelect', 
                      condition = rv.track$typeSelect == 'Random')
      shinyjs::toggle(id = 'colSelect', 
                      condition = rv.track$typeSelect == 'Column')

    })


    observeEvent(input$rst_btn, {
      rv.track$typeSelect <- "None"
      rv.track$listSelect <- NULL
      rv.track$randSelect <- ''
      rv.track$colSelect <- NULL
      rv.track$indices <- NULL
    })

    # Catch event on the list selection
    observeEvent(input$listSelect,{
      rv.track$listSelect <- input$listSelect
      rv.track$randSelect <- ''
      rv.track$colSelect <- ''

      if(is.null(rv.track$listSelect))
        rv.track$indices <- NULL
      else
        rv.track$indices <-  match(rv.track$listSelect, rowData(object())[[idcol(object())]])
    })





    observeEvent(input$randSelect, {
      rv.track$randSelect <- input$randSelect
      rv.track$listSelect  <- NULL
      rv.track$colSelect <- NULL
      cond <- is.null(rv.track$randSelect) 
      cond <- cond || rv.track$randSelect == '' 
      cond <- cond || (as.numeric(rv.track$randSelect) < 0)
      if (!cond)
        rv.track$indices <- sample(seq_len(nrow(object())), 
                                   as.numeric(rv.track$randSelect), 
                                   replace = FALSE)

    })

    observeEvent(input$colSelect, {
      rv.track$colSelect <- input$colSelect
      rv.track$listSelect <- NULL
      rv.track$randSelect <- ''

      if (rv.track$colSelect != '')
        rv.track$indices <- which(rowData(object())[,rv.track$colSelect] == 1)
     })

    return(reactive({rv.track$indices}))

  })

}

## To be copied in the UI
# mod_plots_tracking_ui("plots_tracking_ui_1")

## To be copied in the server
# callModule(mod_plots_tracking_server, "plots_tracking_ui_1")

