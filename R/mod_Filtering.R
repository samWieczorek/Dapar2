
#' @title Module filtering
#' 
#' @param id
#' 
#' @rdname mod_Filtering
#' 
#' @export
#' 
mod_Filtering_ui <- function(id){
  ns <- NS(id)
}


#' @param id xxx
#' @param dataIn xxx
#' @param steps.enabled xxx
#' @param remoteReset xxx
#' @param steps.status xxx
#' @param current.pos xxx
#' @param verbose xxx
#' 
#' @rdname mod_Filtering
#' @importFrom shinyjs toggle hidden
#' @importFrom DT dataTableOutput
#' @importFrom stats setNames
#' @importFrom Magellan Timestamp toggleWidget Get_Worflow_Core_Code
#' 
#' @export
#' 
mod_Filtering_server <- function(id,
                                  dataIn = reactive({NULL}),
                                  steps.enabled = reactive({NULL}),
                                  remoteReset = reactive({FALSE}),
                                  steps.status = reactive({NULL}),
                                  current.pos = reactive({1}),
                                  verbose = FALSE
){
  
  
  # This list contains the basic configuration of the process
  config <- list(
    # Define the type of module
    mode = 'process',
    
    # List of all steps of the process
    steps = c('Description', 'Quanti metadata filtering', 'String based filtering', 'Numerical filtering', 'Save'),
    
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, FALSE, FALSE, FALSE, TRUE),
    
    path_to_md_dir =  system.file('md/', package='DaparToolshed')
  )
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    MetacellTag = "None",
    MetacellFilters = "None",
    KeepRemove = 'delete',
    metacell_value_th = 0,
    choose_metacell_percent_th = 0,
    metacell_value_percent = 0,
    val_vs_percent = 'Value',
    metacellFilter_operator = '<='

  )
  
  
  
  
  
  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(str2expression(Get_Worflow_Core_Code(
      w.names = names(widgets.default.values)
    )))
    
    rv.custom <- reactiveValues(
      temp.filtered = NULL,
      deleted.stringBased = NULL,
      deleted.metacell = NULL,
      deleted.numeric = NULL,
      DT_filterSummary = data.frame(Filter=NULL,
                                     Prefix=NULL,
                                     nbDeleted=NULL,
                                     Total=NULL,
                                     stringsAsFactors=F),
      DT_numfilterSummary = data.frame(Filter=NULL,
                 Condition=NULL,
                 nbDeleted=NULL,
                 Total=NULL,
                 stringsAsFactors=F),
      metacell_Filter_SummaryDT = data.frame(query = NULL,
                                              nbDeleted=NULL,
                                              Total=NULL,
                                              stringsAsFactors=F)
    )
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      file <- paste0(config$path_to_md_dir, '/', id, '.md')
      
      tagList(
        # In this example, the md file is found in the module_examples directory
        # but with a real app, it should be provided by the package which
        # contains the UI for the different steps of the process module.
        # system.file(xxx)
        
        if (file.exists(file))
          includeMarkdown(file)
        else
          p('No Description available'),
        
        
        # Used to show some information about the dataset which is loaded
        # This function must be provided by the package of the process module
        uiOutput(ns('datasetDescription_ui')),
        
        # Insert validation button
        uiOutput(ns('Description_btn_validate_ui'))
      )
    })
    
    output$datasetDescription_ui <- renderUI({
      # Insert your own code to vizualise some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"),
                             "Start",
                             class = btn_success_color)
      Magellan::toggleWidget(rv$steps.enabled['Description'], widget)
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      rv.custom$temp.filtered <- dataIn()
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- global$VALIDATED
    })
    
    
    
    #-------------------------------------------------------------------
    output$Quantimetadatafiltering <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by Magellan
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        div(
          # id = "screen1Filtering",
          
          mod_query_metacell_ui(ns('query')),
          tags$hr(),
          div( style="display:inline-block; vertical-align: middle; align: center;",
               DT::dataTableOutput(ns("metacell_Filter_SummaryDT"))
          ),
          
          hr(),
          ################## Plots section #############################
          mod_plotsMetacellHistos_ui(ns("MVPlots_filtering"))
        ),
        
        # Insert validation button
        uiOutput(ns('Quantimetadatafiltering_btn_validate_ui'))
      )
    })
    
    
    mod_plotsMetacellHistos_server(id = "MVPlots_filtering", 
                                   obj = reactive({rv$current.obj}),
                                   pal = reactive({rv$PlotParams$paletteForConditions}),
                                   pattern = reactive({rv$widgets$filtering$MetacellTag})
    )
    
    
    indices <- mod_query_metacell_server(id = 'query',
                                         obj = reactive({rv$dataIn}),
                                         list_tags = reactive({c('None' = 'None',
                                                                 qMetadata.def(typeDataset(rv$current.obj))$node
                                         )}),
                                         keep_vs_remove = reactive({setNames(nm = c("delete", "keep"))}),
                                         filters = reactive({c("None" = "None",
                                                               "Whole Line" = "WholeLine",
                                                               "Whole matrix" = "WholeMatrix",
                                                               "For every condition" = "AllCond",
                                                               "At least one condition" = "AtLeastOneCond")}),
                                         val_vs_percent = reactive({setNames(nm=c('Count', 'Percentage'))}),
                                         operator = reactive({setNames(nm = SymFilteringOperators())})
    )
    
    
    
    
    
    
    output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
                             "Perform metacell filtering",
                             class = btn_success_color)
      Magellan::toggleWidget(rv$steps.enabled['Quantimetadatafiltering'], widget)
    })
    
    
    observeEvent(input$Quantimetadatafiltering_btn_validate, {
      print('Perform')
      # rv$dataIn <- RunAggregation()
      
      
      
      # rv.custom$temp.agregate <- RunAggregation()
      # 
      # if(is.null(rv.custom$temp.agregate$issues)){
      #   dataOut$trigger <- Magellan::Timestamp()
      #   dataOut$value <- rv$dataIn
      #   rv$steps.status['Agregation'] <- global$VALIDATED
      # } else {
      #   print("error")
      #   shinyjs::toggle('downloadAggregationIssues', 
      #                   condition = !rvModProcess$moduleAggregationDone[1] && length(rv.custom$temp.agregate$issues) > 0
      #   )
      # }
      
      
    })
    

    
    #-------------------------------------------------------
    
    output$StringBasedFiltering <- renderUI({
      wellPanel(
        uiOutput(ns('StringBasedFiltering_btn_validate_ui'))
      )
    })
    
    
    
    
    output$StringBasedFiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("StringBasedFiltering_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      Magellan::toggleWidget(rv$steps.enabled['StringBasedFiltering'], widget)
    })
    
    
    observeEvent(input$StringBasedFiltering_btn_validate, {
     
      
      #dataOut$trigger <- Magellan::Timestamp()
      #dataOut$value <- rv$dataIn
      rv$steps.status['AddMetadata'] <- global$VALIDATED
    })
    
    
    
    
    
    #--------------------------------------------------------
    output$NumericalFiltering <- renderUI({
      wellPanel(
        uiOutput(ns('NumericalFiltering_btn_validate_ui'))
      )
    })
    
    output$NumericalFiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("NumericalFiltering_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      Magellan::toggleWidget(rv$steps.enabled['NumericalFiltering'], widget)
    })
    
    
    observeEvent(input$NumericalFiltering_btn_validate, {
      
      
      
      #dataOut$trigger <- Magellan::Timestamp()
      #dataOut$value <- rv$dataIn
      rv$steps.status['NumericalFiltering'] <- global$VALIDATED
    })
    
    
    
    
    
    #--------------------------------------------------------
    output$Save <- renderUI({
      wellPanel(
        uiOutput(ns('Save_btn_validate_ui'))
      )
    })
    
    output$Save_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Save_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      Magellan::toggleWidget(rv$steps.enabled['Save'], widget)
    })
    
    
    observeEvent(input$Save_btn_validate, {
      
      
      
      #dataOut$trigger <- Magellan::Timestamp()
      #dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
    })
    
    
    #----------------------------------------------------
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

