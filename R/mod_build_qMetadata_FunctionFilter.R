#' @title xxxx
#' 
#' @description xxx
#' 
#' @name query_qMetadata
#' 
#' @value A `list()` of three items:
#' - a trigger : xxx
#' - an instance of the class `FunctionFilter`
#' - A `character()` which describes the query in natural language
#'
#' @examples 
#' if (interactive()){
#' data(ft_na)
#' ui <- mod_build_qMetadata_FunctionFilter_ui('query')
#' 
#' server <- function(input, output, session) {
#'   
#'   rv <- reactiveValues(
#'     res = NULL
#'   )
#'   ll.tags <- c('None' = 'None', 
#'                qMetadata.def(typeDataset(ft_na[[1]]))$node)
#'   rv$res <- mod_build_qMetadata_FunctionFilter_server('query', 
#'                                                       obj = reactive({ft_na[[1]]}),
#'                                                       conds = reactive({colData(ft_na)$Condition}),
#'                                                       list_tags = reactive({ll.tags}),
#'                                                       keep_vs_remove = reactive({setNames(nm = c("delete", "keep"))}),
#'                                                       val_vs_percent = reactive({setNames(nm=c('Count', 'Percentage'))}),
#'                                                       operator = reactive({setNames(nm = SymFilteringOperators())})
#'   )
#'   
#'   
#'   observeEvent(rv$res$dataOut()$trigger, ignoreNULL = TRUE, ignoreInit = TRUE, {
#'     print(rv$res$dataOut()$fun)
#'   })
#' }
#' 
#' shinyApp(ui=ui, server=server)
#' }
NULL









#' @param id xxx
#' 
#' @export
#' 
#' @rdname query_qMetadata
#' 
mod_build_qMetadata_FunctionFilter_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    div(
      fluidRow(
        column(2, uiOutput(ns('chooseTag_ui'))),
        column(2, uiOutput(ns("chooseKeepRemove_ui"))),
        column(2, uiOutput(ns('chooseScope_ui'))),
        column(6, tagList(
          uiOutput(ns("qMetadataScope_widgets_set2_ui")))
        )
      ),
      div( style="display:inline-block; vertical-align: middle; align: center;",
           uiOutput(ns('qMetadataScope_request_ui'))
      ),
      actionButton(ns("BuildFilter_btn"), "Build filter")
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
mod_build_qMetadata_FunctionFilter_server <- function(id,
                                      obj,
                                      conds,
                                      list_tags = reactive({NULL}),
                                      keep_vs_remove = reactive({NULL}),
                                      val_vs_percent = reactive({NULL}),
                                      operator = reactive({NULL}),
                                      reset = reactive({NULL}),
                                      is.enabled = reactive({TRUE})
                                      ) {
  
  rv <- reactiveValues(
    indices = NULL,
    trigger = NULL,
    functionFilter = NULL,
    query = NULL
  )
  
  moduleServer(id,function(input, output, session) {
    ns <- session$ns
    mod_helpPopover_server("tag_help", 
                          title = "Nature of data to filter", 
                          content = "Define xxx")
    help.txt1 <- "To filter the missing values, the choice of the lines to be kept is 
    made by different options:
    <ul>
    <li><strong>None</strong>: No filtering, the quantitative data is left unchanged.</li>
    <li><strong>(Remove) Empty lines</strong>: All the lines with 100% of missing 
    values are filtered out.</li>
    <li><strong>Whole Matrix</strong>: The lines (across all conditions) which contain 
    less quantitative value than a user-defined threshold are kept;</li>
    <li><strong>For every condition</strong>: The lines for which each condition contain 
    less quantitative value than a user-defined threshold are deleted;</li>
    <li><strong>At least one condition</strong>: The lines for which at least one condition 
    contain less quantitative value than a user-defined threshold are deleted.</li>
    </ul>"
      
    mod_helpPopover_server("filterScope_help", 
                          title = "Scope", 
                          content = HTML(help.txt1)
                          )
    
    rv.widgets <- reactiveValues(
      tag = "None",
      scope = "None",
      keepRemove = 'delete',
      valueTh = 0,
      percentTh = 0,
      valPercent = 'Count',
      operator = '<='
      )
    
    dataOut <- reactiveValues(
      trigger = NULL,
      fun = NULL,
      query = NULL
    )
    
    observeEvent(reset(),{
      rv.widgets$tag <- "None"
      rv.widgets$scope <- "None"
      rv.widgets$keepRemove <- 'delete'
      rv.widgets$valueTh <- 0
      rv.widgets$percentTh <- 0
      rv.widgets$valPercent <- 'Count'
      rv.widgets$operator <- '<='
      updateSelectInput(session, "chooseTag", selected = "None")
      })
    
    observeEvent(input$chooseTag, { rv.widgets$tag <- input$chooseTag})
    observeEvent(input$chooseKeepRemove, {  rv.widgets$keepRemove <- input$chooseKeepRemove})
    observeEvent(input$chooseScope, {  rv.widgets$scope <- input$chooseScope})
    observeEvent(input$chooseValPercent, {rv.widgets$valPercent <- input$chooseValPercent})
    observeEvent(input$chooseValueTh, { rv.widgets$valueTh <- input$chooseValueTh})
    observeEvent(input$choosePercentTh, {rv.widgets$percentTh <- input$choosePercentTh})
    observeEvent(input$chooseOperator, {  rv.widgets$operator <- input$chooseOperator})
                 
    
    output$chooseTag_ui <- renderUI({
      widget <- selectInput(ns("chooseTag"),
                            mod_helpPopover_ui(ns("tag_help")),
                            choices = list_tags(),
                            selected = rv.widgets$tag,
                            width='200px')
      Magellan::toggleWidget(widget, is.enabled())
      })
    
    output$chooseKeepRemove_ui <- renderUI({
      req(rv.widgets$tag != 'None')
      widget <- radioButtons(ns("chooseKeepRemove"),
                             "Type of filter operation",
                             choices = keep_vs_remove(),
                             selected = rv.widgets$keepRemove)
      Magellan::toggleWidget(widget, is.enabled())
      })
    
    output$chooseScope_ui <- renderUI({
      req(rv.widgets$tag != 'None')
      widget <- selectInput(ns("chooseScope"),
                            mod_helpPopover_ui(ns("filterScope_help")),
                            choices = c("None" = "None",
                                        "Whole Line" = "WholeLine",
                                        "Whole matrix" = "WholeMatrix",
                                        "For every condition" = "AllCond",
                                        "At least one condition" = "AtLeastOneCond"),
                            selected = rv.widgets$scope,
                            width='200px')
      Magellan::toggleWidget(widget, is.enabled())
      })
    
    
    
    
    output$qMetadataScope_widgets_set2_ui <- renderUI({
      req(!(rv.widgets$scope %in% c("None", "WholeLine")))
      mod_helpPopover_server("chooseValPercent_help", 
                             title = paste("#/% of values to ", rv.widgets$keepRemove),
                             content = "Define xxx")
      
      tagList(
        fluidRow(
          column(4,
                 radioButtons(ns('chooseValPercent'),
                              mod_helpPopover_ui(ns("chooseValPercent_help")),
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
                 uiOutput(ns('chooseValue_ui')),
                 uiOutput(ns('choosePercentage_ui'))
                 )
          )
        )
      })
    
    
    output$chooseValue_ui <- renderUI({
      req(rv.widgets$valPercent == 'Count')
      req(!(rv.widgets$scope %in% c("None", "WholeLine")))
      mod_helpPopover_server("value_th_help", 
                             title = "Count threshold", 
                             content = "Define xxx")
      
      tagList(
        mod_helpPopover_ui(ns("keepVal_help")),
        selectInput(ns("chooseValueTh"),
                    mod_helpPopover_ui(ns("value_th_help")),
                    choices = getListNbValuesInLines(object = obj(), 
                                                     conds = conds(), 
                                                     type = rv.widgets$scope
                                                     ),
                    selected = rv.widgets$valueTh,
                    width='150px')
        )
      })
    
    output$choosePercentage_ui <- renderUI({
      req(rv.widgets$valPercent == 'Percentage')
      req(!(rv.widgets$scope %in% c("None", "WholeLine")))
      
      mod_helpPopover_server("percentTh_help", 
                             title = "Percentage threshold", 
                             content = "Define xxx")
      
      tagList(
        mod_helpPopover_ui(ns("keepVal_percent_help")),
        sliderInput(ns("choosePercentTh"), 
                    mod_helpPopover_ui(ns("percentTh_help")),
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
      if (rv.widgets$scope == "None"){
        txt_summary <- "No filtering is processed."
        } else if (rv.widgets$scope == "WholeLine") {
          txt_summary <- paste(rv.widgets$keepRemove,
                               "lines that contain only",
                               rv.widgets$tag)
          } else {
            text_method <-  switch(rv.widgets$scope,
                                   WholeMatrix = "the whole matrix.",
                                   AllCond = "each condition.",
                                   AtLeastOneCond = "at least one condition.")
            
            if(rv.widgets$valPercent == 'Count'){
              text_threshold <- rv.widgets$valueTh
              } else {
                text_threshold <- paste(as.character(rv.widgets$percentTh)," %", sep="")
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
      tags$p(paste("You are going to ", WriteQuery()), 
             style = "font-size: small; text-align : center; color: purple;")
      })
    
    
    #Set useless widgets to default values
    observeEvent(rv.widgets$scope == 'WholeLine',{
      rv.widgets$percentThh <- 0
      rv.widgets$valueTh <- 0
      rv.widgets$valPercent <- 'Percentage'
      },
      priority = 1000)
    
   
    
    
    BuildFunctionFilter <- reactive({
      req(obj())
      req(rv.widgets$tag != 'None')
      req(rv.widgets$scope != 'None')
      
      th <- switch(rv.widgets$valPercent,
                   Percentage = rv.widgets$percentTh / 100,
                   Count = as.integer(rv.widgets$valueTh)
                   )
      
      ff <- switch(rv.widgets$scope,
                   WholeLine = FunctionFilter('qMetadataWholeLine', 
                                              cmd = rv.widgets$keepRemove, 
                                              pattern = rv.widgets$tag),
                   WholeMatrix = FunctionFilter('qMetadataWholeMatrix',
                                                cmd = rv.widgets$keepRemove,
                                                pattern = rv.widgets$tag,
                                                percent = rv.widgets$valPercent,
                                                th = th,
                                                operator = rv.widgets$operator),
                   AllCond = FunctionFilter('qMetadataOnConditions', 
                                            cmd = rv.widgets$keepRemove, 
                                            mode = rv.widgets$scope,
                                            pattern = rv.widgets$tag,
                                            conds = conds(),
                                            percent = rv.widgets$valPercent, 
                                            th = th, 
                                            operator = rv.widgets$operator),
                   AtLeastOneCond = FunctionFilter('qMetadataOnConditions', 
                                                   cmd = rv.widgets$keepRemove, 
                                                   mode = rv.widgets$scope,
                                                   pattern = rv.widgets$tag,
                                                   conds = conds(),
                                                   percent = rv.widgets$valPercent,
                                                   th = th, 
                                                   operator = rv.widgets$operator)
      )

      })
    
    
    observeEvent(input$BuildFilter_btn, {
      dataOut$trigger <- as.numeric(Sys.time())
      dataOut$fun <- BuildFunctionFilter()
      dataOut$query <- WriteQuery()
    })
    
    
    
    list(dataOut = reactive({dataOut}))

    })
}

