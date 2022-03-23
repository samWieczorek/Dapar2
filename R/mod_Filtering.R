
#' @title Module Filtering
#' 
#' @description 
#' 
#' @section Step 'Quanti metadata filtering':
#' 
#' xxxxxxx
#' 
#' @section Step 'Variable filtering':
#' 
#' xxxxx
#' 
#' @section Step 'Save':
#' 
#' xxxxxx
#' 
#' @seealso The user manual of the package `Magellan`.
#'
#' 
#' @name mod_filtering
#' 
#' @examples 
#' if(interactive()){
#'  run_workflow('Filtering', verbose = TRUE)
#' }
NULL



#' @param id
#' 
#' @rdname mod_Filtering
#' 
#' @export
#' 
mod_Filtering_ui <- function(id){
  ns <- NS(id)
}


#' @param id A `character(1)` which is the 'id' of the module.
#' @param dataIn An instance of the class `QFeatures`
#' @param steps.enabled A `logical()` which indicates whether each step is
#' enabled or disabled in the UI.
#' @param remoteReset A `logical(1)` which acts asa a remote command to reset
#' the module to its default values. Default is FALSE.
#' @param steps.status A `logical()` which indicates the status of each step
#' which can be either 'validated', 'undone' or 'skipped'.
#' enabled or disabled in the UI.
#' @param current.pos A `interger(1)` which acts as a remote command to make
#'  a step active in the timeline. Default is 1.
#' @param verbose A `logical(1)` to indicates whether run and debug infos must
#' be printed in the console. Default is FALSE.
#' 
#' @rdname mod_Filtering
#' @importFrom shinyjs toggle hidden
#' @importFrom DT dataTableOutput renderDataTable DTOutput datatable
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
  config <- Magellan::Config(
    name = 'Filtering',
    # Define the type of module
    mode = 'process',
    
    # List of all steps of the process
    steps = c('Description', 
              'Quanti metadata filtering', 
              'Variable filtering', 
              'Save'),
    
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, FALSE, FALSE, TRUE),
    
    path_to_md_dir =  system.file('md/', package='DaparToolshed')
  )
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    tag = "None",
    scope = "None",
    keepRemove = 'delete',
    valueTh = 0,
    percentTh = 0,
    valuePercent = 0,
    valPercent = 'Value',
    operator = '<='
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
      funFilter = NULL,
      
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
      qMetadata_Filter_SummaryDT = data.frame(query = NULL,
                                              nbDeleted = NULL,
                                              TotalMainAssay = NULL,
                                              stringsAsFactors = F)
    )
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      file <- paste0(config@path_to_md_dir, '/', id, '.md')
      
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
      Magellan::toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- global$VALIDATED
    })
    
    
    
    #--------------------------------------------------------------
    # Quantitative metadata filtering UI
    #--------------------------------------------------------------
    
    output$Quantimetadatafiltering <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by Magellan
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
       # mod_build_qMetadata_FunctionFilter_ui(ns('query')),
        DT::dataTableOutput(ns("qMetadata_Filter_Summary_DT")),
        uiOutput(ns('Quantimetadatafiltering_buildQuery_ui')),
        
        uiOutput(ns('example_ui')),
        mod_ds_qMetadata_ui(ns('plots')),
        # Insert validation button
        uiOutput(ns('Quantimetadatafiltering_btn_validate_ui'))
      )
    })
    
    
    
    output$example_ui <- renderUI({
      req(rv.custom$funFilter$dataOut()$fun)
      req(rv$steps.status['Quantimetadatafiltering'] == 0)
      
      temp <- filterFeaturesOneSE(object = rv$dataIn[[length(rv$dataIn)]],
                                  filters = list(rv.custom$funFilter$dataOut()$fun)
                                  )
      
      mod_filtering_example_server(id = 'filteringExample',
                                   objBefore = reactive({mainAssay(rv$dataIn)}),
                                   objAfter = reactive({temp}),
                                   query = reactive({rv.custom$funFilter$dataOut()$query})
                                   )
      widget <- mod_filtering_example_ui(ns('filteringExample'))
      Magellan::toggleWidget(widget, rv$steps.enabled['Quantimetadatafiltering'])

    })
    
    
    mod_ds_qMetadata_server(id = 'plots',
                            se = reactive({mainAssay(rv$dataIn)}),
                            init.pattern = 'missing',
                            conds = design(rv$dataIn)$Condition
                            )
      
    
    output$qMetadata_Filter_Summary_DT <- DT::renderDataTable(server=TRUE,{
      req(rv$dataIn)
      req(rv.custom$qMetadata_Filter_SummaryDT)
      isolate({
        
        if (nrow(rv.custom$qMetadata_Filter_SummaryDT )==0){
          df <- data.frame(query = "-",
                           nbDeleted = 0,
                           TotalMainAssay = nrow(mainAssay(rv$dataIn)),
                           stringsAsFactors = FALSE)
          rv.custom$qMetadata_Filter_SummaryDT <- df
        }
        
        
        DT::datatable(rv.custom$qMetadata_Filter_SummaryDT,
                      extensions = c('Scroller'),
                      rownames = FALSE,
                      options = list(
                        dom = 'rt',
                        initComplete = .initComplete(),
                        deferRender = TRUE,
                        bLengthChange = FALSE
                      ))
      })
    })
    
    
    
    
    
    
    output$Quantimetadatafiltering_buildQuery_ui <- renderUI({
      widget <- mod_build_qMetadata_FunctionFilter_ui(ns("query"))
      Magellan::toggleWidget(widget, rv$steps.enabled['Quantimetadatafiltering'])
    })
    
    rv.custom$funFilter <- mod_build_qMetadata_FunctionFilter_server(id = 'query',
                                         obj = reactive({mainAssay(rv$dataIn)}),
                                         conds = reactive({design(rv$dataIn)$Condition}),
                                         list_tags = reactive({
                                           req(rv$dataIn)
                                           c('None' = 'None',
                                             qMetadata.def(typeDataset(mainAssay(rv$dataIn)))$node)}),
                                         keep_vs_remove = reactive({setNames(nm = c("delete", "keep"))}),
                                         val_vs_percent = reactive({setNames(nm=c('Count', 'Percentage'))}),
                                         operator = reactive({setNames(nm = SymFilteringOperators())})
    )

    
    output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
                             "Perform qMetadata filtering",
                             class = btn_success_color)
      cond <- !is.null(rv.custom$funFilter$dataOut()$fun) && rv$steps.enabled['Quantimetadatafiltering']
      Magellan::toggleWidget(widget, cond)
    })

    observeEvent(input$Quantimetadatafiltering_btn_validate, {
      
      rv$dataIn <- filterFeaturesOneSE(object = rv$dataIn,
                                       i = length(rv$dataIn),
                                       name = 'qMetadataFiltered',
                                       filters = list(rv.custom$funFilter$dataOut()$fun)
      )

      # Add infos
       nBefore <- nrow(rv$dataIn[[length(rv$dataIn) - 1]])
       nAfter <- nrow(rv$dataIn[[length(rv$dataIn)]])
       
      
      df <- data.frame(query =  rv.custom$funFilter$dataOut()$query,
                       nbDeleted = nBefore - nAfter,
                       TotalMainAssay = nrow(rv$dataIn[[length(rv$dataIn)]]))
      
      rv.custom$qMetadata_Filter_SummaryDT <- rbind(rv.custom$qMetadata_Filter_SummaryDT , 
                                                     df)
      
      params(rv$dataIn, length(rv$dataIn)) <- reactiveValuesToList(rv.widgets)
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Quantimetadatafiltering'] <- global$VALIDATED
      
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
      Magellan::toggleWidget(widget, rv$steps.enabled['StringBasedFiltering'])
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
      Magellan::toggleWidget(widget, rv$steps.enabled['NumericalFiltering'])
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
      Magellan::toggleWidget(widget, rv$steps.enabled['Save'])
    })
    
    
    observeEvent(input$Save_btn_validate, {
      #browser()
      # Add the parameters values to the new dataset
      par <- lapply(names(widgets.default.values),
                        function(x) x )
      
      
      params(rv$dataIn[[length(rv$dataIn)]]) <- par
      
    dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
    })
    
    
    #----------------------------------------------------
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

