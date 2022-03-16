#' @title xxxx
#' 
#' @description xxx
#' 
#' @name query_qMetadata
#'
#' @examples 
#' if (interactive()){
#' data(ft_na)
#'  ui <- mod_query_qMetadata_ui('query')
#' 
#'  server <- function(input, output, session) {
#'  ll.tags <- c('None' = 'None', qMetadata.def(typeDataset(ft_na[[1]]))$node)
#'   mod_query_qMetadata_server('query', 
#'                              obj = reactive({ft[[1]]}),
#'                              list_tags = reactive({ll.tags}),
#'                              keep_vs_remove = reactive({setNames(nm = c("delete", "keep"))}),
#'                              filters = reactive({c("None" = "None",
#'                              "Whole Line" = "WholeLine",
#'                              "Whole matrix" = "WholeMatrix",
#'                              "For every condition" = "AllCond",
#'                              "At least one condition" = "AtLeastOneCond")}),
#'                              val_vs_percent = reactive({setNames(nm=c('Count', 'Percentage'))}),
#'                              operator = reactive({setNames(nm = SymFilteringOperators())})
#'                              )
#'                              }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
NULL









#' @param id xxx
#' 
#' @export
#' 
#' @rdname query_qMetadata
#' 
mod_query_qMetadata_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    div(
      fluidRow(
        column(2, uiOutput(ns('chooseTag_ui'))),
        column(2, uiOutput(ns("ChoosekeepRemove_ui"))),
        column(2, uiOutput(ns('chooseFilters_ui'))),
        column(6, tagList(
          uiOutput(ns("example_ui")),
          uiOutput(ns("qMetadataFilters_widgets_set2_ui")))
        )
      ),
      div( style="display:inline-block; vertical-align: middle; align: center;",
           uiOutput(ns('qMetadataFilter_request_ui'))
      )
    )
  )
  
}




#' @param id xxx
#' @param obj xxx
#' @param list_tags xxx
#' @param keep_vs_remove xxx
#' @param filters xxx
#' @param val_vs_percent xxx
#' @param operator xxx
#' @param reset xxx
#' 
#' @name mod_query_qMetadata_server
#' @rdname quantitative-metadata
#' 
#' @export
#' 
mod_query_qMetadata_server <- function(id,
                                      obj,
                                      list_tags = reactive({NULL}),
                                      keep_vs_remove = reactive({NULL}),
                                      filters = reactive({NULL}),
                                      val_vs_percent = reactive({NULL}),
                                      operator = reactive({NULL}),
                                      reset = reactive({NULL})
) {
  
  rv <- reactiveValues(
    indices = NULL,
    trigger = NULL,
    params = NULL,
    query = NULL
  )
  
  moduleServer(id,function(input, output, session) {
    ns <- session$ns
    mod_popover_for_help_server("tag_help", 
                            data = reactive(list(title = "Nature of data to filter", 
                                                 content="Define xxx")))
                 
    mod_popover_for_help_server("filterScope_help", 
                            data = reactive(list(title = "Scope", 
                                                 content=HTML(paste0("To filter the missing values, the choice of the lines to be kept is made by different options:"),
                                                              ("<ul>"),
                                                              ("<li><strong>None</strong>: No filtering, the quantitative data is left unchanged.</li>"),
                                                              ("<li><strong>(Remove) Empty lines</strong>: All the lines with 100% of missing values are filtered out.</li>"),
                                                              ("<li><strong>Whole Matrix</strong>: The lines (across all conditions) which contain less quantitative value than a user-defined threshold are kept;</li>"),
                                                              ("<li><strong>For every condition</strong>: The lines for which each condition contain less quantitative value than a user-defined threshold are deleted;</li>"),
                                                              ("<li><strong>At least one condition</strong>: The lines for which at least one condition contain less quantitative value than a user-defined threshold are deleted.</li>"),
                                                              ("</ul>")
                                                 )
                            )))
    rv.widgets <- reactiveValues(
      tag = "None",
      filters = "None",
      KeepRemove = 'delete',
      value_th = 0,
      percent_th = 0,
      val_vs_percent = 'Count',
      operator = '<='
      )
    
    observeEvent(reset(),{
      rv.widgets$tag <- "None"
      rv.widgets$filters <- "None"
      rv.widgets$keepRemove <- 'delete'
      rv.widgets$valueTh <- 0
      rv.widgets$percentTh <- 0
      rv.widgets$valPercent <- 'Count'
      rv.widgets$operator <- '<='
      updateSelectInput(session, "chooseTag", selected = "None")
      })
    
    observeEvent(input$chooseTag, { rv.widgets$tag <- input$chooseTag})
    observeEvent(input$chooseKeepRemove, {  rv.widgets$keepRemove <- input$chooseKeepRemove})
    observeEvent(input$chooseFilters, {  rv.widgets$filters <- input$chooseFilters})
    observeEvent(input$chooseValPercent, {rv.widgets$valPercent <- input$chooseValPercent})
    observeEvent(input$chooseValueTh, { rv.widgets$alueTh <- input$chooseValueTh})
    observeEvent(input$choosePercentTh, {rv.widgets$percentTh <- input$choosePercentTh})
    observeEvent(input$chooseOperator, {  rv.widgets$operator <- input$chooseOperator})
                 
    
    output$chooseTag_ui <- renderUI({
      selectInput(ns("chooseTag"),
                  mod_popover_for_help_ui(ns("tag_help")),
                  choices = list_tags(),
                  selected = rv.widgets$tag,
                  width='200px')
      })
    
    output$chooseKeepRemove_ui <- renderUI({
      req(rv.widgets$tag != 'None')
      radioButtons(ns("ChooseKeepRemove"),
                   "Type of filter operation",
                   choices = keep_vs_remove(),
                   selected = rv.widgets$keepRemove)
      })
    
    output$chooseFilters_ui <- renderUI({
      req(rv.widgets$tag != 'None')
      selectInput(ns("chooseFilters"),
                  mod_popover_for_help_ui(ns("filterScope_help")),
                  choices = filters(),
                  selected = rv.widgets$filters,
                  width='200px')
      })
    
    output$example_ui <- renderUI({
      req(rv.widgets$filters != "None")
      mod_filtering_example_server(id = 'filteringExample',
                                   obj = reactive({obj()}),
                                   indices = reactive({CompileIndices()}),
                                   params = reactive({rv.widgets}),
                                   txt = reactive({WriteQuery()})
                                   )
      
      mod_filtering_example_ui(ns('filteringExample'))
      })
    
    output$qMetadataFilters_widgets_set2_ui <- renderUI({
      req(!(rv.widgets$filters %in% c("None", "WholeLine")))
      mod_popover_for_help_server("choose_val_vs_percent_help", 
                                  data = reactive(list(title = paste("#/% of values to ", rv.widgets$KeepRemove),
                                                   content="Define xxx")))
      
      tagList(
        fluidRow(
          column(4,
                 radioButtons(ns('chooseValPercent'),
                              mod_popover_for_help_ui(ns("chooseValPercent_help")),
                              choices = val_vs_percent(),
                              selected = rv.widgets$valPercent
                              )
                 ),
          column(8,
                 selectInput(ns("chooseOperator"),
                             "Choose operator",
                             choices = operator(),
                             selected = rv.widgets$operator,
                             width='100px'),
                 uiOutput(ns('choose_value_ui')),
                 uiOutput(ns('choose_percentage_ui'))
                 )
          )
        )
      })
    output$choose_value_ui <- renderUI({
      req(rv.widgets$val_vs_percent == 'Count')
      req(!(rv.widgets$qMetadataFilters %in% c("None", "WholeLine")))
      mod_popover_for_help_server("value_th_help", 
                              data = reactive(list(title = "Count threshold", 
                                                   content="Define xxx")))
      
      tagList(
        mod_popover_for_help_ui(ns("keepVal_help")),
        selectInput(ns("chooseValueTh"),
                    mod_popover_for_help_ui(ns("value_th_help")),
                    choices = getListNbValuesInLines(obj(), 
                                                     type = rv.widgets$filters),
                    selected = rv.widgets$valueTh,
                    width='150px')
        )
      })
    
    output$choosePercentage_ui <- renderUI({
      req(rv.widgets$valPercent == 'Percentage')
      req(!(rv.widgets$filters %in% c("None", "WholeLine")))
      
      mod_popover_for_help_server("percentTh_help", 
                              data = reactive(list(title = "Percentage threshold", 
                                                   content="Define xxx")))
      
      tagList(
        mod_popover_for_help_ui(ns("keepVal_percent_help")),
        sliderInput(ns("choosePercentTh"), 
                    mod_popover_for_help_ui(ns("percentTh_help")),
                    min = 0,
                    max = 100,
                    step = 1,
                    value = rv.widgets$percentTh,
                    width='250px')
        )
      })
    
    
    ## ---------------------------------------------------------------------
    ##
    ## This function xxx
    ##
    ## ---------------------------------------------------------------------
    WriteQuery <- reactive({
      if (rv.widgets$filters == "None"){
        txt_summary <- "No filtering is processed."
        } else if (rv.widgets$filters == "WholeLine") {
          txt_summary <- paste(rv.widgets$keepRemove,
                               "lines that contain only",
                               rv.widgets$tag)
          } else {
            text_method <-  switch(rv.widgets$filters,
                                   WholeMatrix = "the whole matrix.",
                                   AllCond = "each condition.",
                                   AtLeastOneCond = "at least one condition.")
            
            if(rv.widgets$valPercent == 'Count'){
              text_threshold <- rv.widgets$valueTh
              } else {
                text_threshold <- paste(as.character(rv.widgets$percentTh),
                                               " %", sep="")
                }
            
            txt_summary <- paste(rv.widgets$keepRemove,
                                 " lines where number of ",
                                 rv.widgets$tag,
                                 " data ",
                                 rv.widgets$operator,
                                 " ",
                                 text_threshold,
                                 " in ",
                                 text_method)
            }
      txt_summary
      })
    
    
    output$qMetadataFilter_request_ui <- renderUI({
      txt_summary <- paste("You are going to ", WriteQuery())
      tags$p(txt_summary, style = "font-size: small; text-align : center; color: purple;")
      })
    
    
    #Set useless widgets to default values
    observeEvent(rv.widgets$filters == 'WholeLine',{
      rv.widgets$percentThh <- 0
      rv.widgets$valueTh <- 0
      rv.widgets$valPercent <- 'Percentage'
      },
      priority = 1000)
    
    
    CompileIndices <- reactive({
      req(obj())
      req(rv.widgets$tag != 'None')
      req(rv.widgets$filters != 'None')
      
      th <- switch(rv.widgets$valPercent,
                   Percentage =  rv.widgets$percentTh / 100,
                   Count = as.integer(rv.widgets$valueTh)
                   )
      
      GetIndices_filtering(obj = obj(),
                           level = typeDataset(obj()),
                           pattern = rv.widgets$tag,
                           type = rv.widgets$filters,
                           percent = rv.widgets$valPercent == 'Percentage',
                           op = rv.widgets$operator,
                           th = th)
      })
    
    reactive({list(trigger = as.numeric(Sys.time()),
                   indices = CompileIndices(),
                   params = list(tag = rv.widgets$tag,
                                 keepRemove = rv.widgets$keepRemove,
                                 filters = rv.widgets$filters,
                                 percentTh = rv.widgets$percentTh,
                                 valueTh = rv.widgets$valueTh,
                                 valPercent = rv.widgets$valPercent,
                                 operator = rv.widgets$operator),
                   query = WriteQuery()
                   )
      })
    })
}

