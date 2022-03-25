
#' @title Module Filtering
#' 
#' @description 
#' 
#' xxxxx
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
    Quantimetadatafiltering_tag = "None",
    Quantimetadatafiltering_scope = "None",
    Quantimetadatafiltering_keepRemove = 'delete',
    Quantimetadatafiltering_valueTh = 0,
    Quantimetadatafiltering_percentTh = 0,
    Quantimetadatafiltering_valuePercent = 0,
    Quantimetadatafiltering_valPercent = 'Value',
    Quantimetadatafiltering_operator = '<=',
    Variablefiltering_cname = 'None',
    Variablefiltering_value = NULL,
    Variablefiltering_operator = ''
  )
  
  
  rv.custom.default.values <- list(
    funFilter = NULL,
    varFilters = list(),
    Variablefiltering_query = list()
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
      w.names = names(widgets.default.values),
      rv.custom.names = names(rv.custom.default.values)
    )))
    
    
    
    
    
    # Add observer to catch the remoteReset so as to 
    # manage rv.custom values
    
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
    #
    #             Quantitative metadata filtering UI
    # 
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
      req(rv.custom$funFilter()$fun)
      req(rv$steps.status['Quantimetadatafiltering'] == 0)
      
      temp <- filterFeaturesOneSE(object = rv$dataIn[[length(rv$dataIn)]],
                                  filters = list(rv.custom$funFilter()$fun)
                                  )
      
      mod_filtering_example_server(id = 'filteringExample',
                                   objBefore = reactive({mainAssay(rv$dataIn)}),
                                   objAfter = reactive({temp}),
                                   query = reactive({rv.custom$funFilter()$query})
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

    
    # Update the list of queries each time a new filter is added
    observeEvent(rv.custom$funFilter()$trigger, {
      
      print(rv.custom$funFilter()$value()$ll.query)
    })
    
    output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
                             "Perform qMetadata filtering",
                             class = btn_success_color)
      cond <- length(rv.custom$funFilter()$value()$ll.fun)
      cond <- cond && rv$steps.enabled['Quantimetadatafiltering']
      Magellan::toggleWidget(widget, cond)
    })

    
    observeEvent(input$Quantimetadatafiltering_btn_validate, {
      
      rv$dataIn <- filterFeaturesOneSE(object = rv$dataIn,
                                       i = length(rv$dataIn),
                                       name = 'qMetadataFiltered',
                                       filters = list(rv.custom$funFilter()$value()$ll.fun)
      )

      # Add infos
       nBefore <- nrow(rv$dataIn[[length(rv$dataIn) - 1]])
       nAfter <- nrow(rv$dataIn[[length(rv$dataIn)]])
       
      
      df <- data.frame(query =  rv.custom$funFilter()$value()$ll.query,
                       nbDeleted = nBefore - nAfter,
                       TotalMainAssay = nrow(rv$dataIn[[length(rv$dataIn)]]))
      
      rv.custom$qMetadata_Filter_SummaryDT <- rbind(rv.custom$qMetadata_Filter_SummaryDT , 
                                                     df)
      par <- rv.custom$funFilter()$value()$ll.widgets.value
      params(rv$dataIn, length(rv$dataIn)) <- par
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Quantimetadatafiltering'] <- global$VALIDATED
      
    })
    

    
    # -------------------------------------------------------
    #
    #                  Variable filtering
    #
    # -------------------------------------------------------
    output$Variablefiltering <- renderUI({
      wellPanel(
        uiOutput(ns('Variablefiltering_cname_ui')),
        uiOutput(ns('Variablefiltering_value_ui')),
        uiOutput(ns('Variablefiltering_operator_ui')),
        
        uiOutput(ns('Variablefiltering_addFilter_btn_ui')),
        
        htmlOutput(ns('show_filters_ui')),
        uiOutput(ns('Variablefiltering_btn_validate_ui')),
        
        DT::dataTableOutput(ns("Variable_Filter_Summary_DT"))
      )

    })
    
    
    
    
    output$show_filters_ui <- renderText({
      req(length(rv.custom$Variablefiltering_query) > 0)
     # browser()
      HTML(paste0(unlist(rv.custom$Variablefiltering_query),
                  collapse = " "))
    })
    
    
    
    output$Variable_Filter_Summary_DT <- DT::renderDataTable(server=TRUE,{
      req(rv$dataIn)
      req(rv.custom$variable_Filter_SummaryDT)
      isolate({
        
        if (nrow(rv.custom$variable_Filter_SummaryDT )==0){
          df <- data.frame(query = "-",
                           nbDeleted = 0,
                           TotalMainAssay = nrow(mainAssay(rv$dataIn)),
                           stringsAsFactors = FALSE)
          rv.custom$variable_Filter_SummaryDT <- df
        }
        
        
        DT::datatable(rv.custom$variable_Filter_SummaryDT,
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
  
  output$Variablefiltering_cname_ui <- renderUI({
    .choices <- c("None", 
                  colnames(SummarizedExperiment::rowData(mainAssay(rv$dataIn))))
  
    widget <- selectInput(ns("Variablefiltering_cname"), 
                          "Column name", 
                          choices = setNames(.choices, nm = .choices),
                          width = '300px')
    
    Magellan::toggleWidget(widget, rv$steps.enabled['Variablefiltering'])
  })
  
  
  output$Variablefiltering_operator_ui <- renderUI({
    req(rv.widgets$Variablefiltering_value)
    if (is.na(as.numeric(rv.widgets$Variablefiltering_value)))
      .operator <- c('==', '!=', 'startsWith', 'endsWith', 'contains')
    else
      .operator <- DaparToolshed::SymFilteringOperators()
    
    
  widget <- selectInput(ns("Variablefiltering_operator"), 
              "operator",
              choices = setNames(nm = .operator),
              width='100px')
  Magellan::toggleWidget(widget, rv$steps.enabled['Variablefiltering'])
  })
  
  output$Variablefiltering_value_ui <- renderUI({
    widget <- textInput(ns("Variablefiltering_value"), 
                          "value",
                          width='100px')
    Magellan::toggleWidget(widget, rv$steps.enabled['Variablefiltering'])
  })
  
  
    output$Variablefiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Variablefiltering_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      Magellan::toggleWidget(widget, 
                             rv$steps.enabled['Variablefiltering'] && length(rv.custom$varFilters) > 0)
    })
    
    
    
    
    
    
    output$Variablefiltering_addFilter_btn_ui <- renderUI({
      widget <- actionButton(ns('Variablefiltering_addFilter_btn'), 'Add filter')
      Magellan::toggleWidget(widget, rv$steps.enabled['Variablefiltering'])
    })
    
    
    observeEvent(input$Variablefiltering_addFilter_btn, {
      
      if (!is.na(as.numeric(rv.widgets$Variablefiltering_value)))
        value <- as.numeric(rv.widgets$Variablefiltering_value)
      else
        rv.widgets$Variablefiltering_value
        
        
      rv.custom$varFilters <- append(rv.custom$varFilters,
                                     VariableFilter(field = rv.widgets$Variablefiltering_cname, 
                                                    value = value,
                                                    condition = rv.widgets$Variablefiltering_operator)
                                     )
      rv.custom$Variablefiltering_query <- append(rv.custom$Variablefiltering_query,
                                                  paste0(rv.widgets$Variablefiltering_cname, " ",
                                                         rv.widgets$Variablefiltering_operator, " ",
                                                         value))
      
    })
    
    
    
    observeEvent(input$Variablefiltering_btn_validate, {
     
      rv$dataIn <- filterFeaturesOneSE(object = rv$dataIn,
                                       i = length(rv$dataIn),
                                       name = 'variableFiltered',
                                       filters = rv.custom$varFilters)
      # Add infos
      nBefore <- nrow(rv$dataIn[[length(rv$dataIn) - 1]])
      nAfter <- nrow(rv$dataIn[[length(rv$dataIn)]])
      
      
      df <- data.frame(query = paste0(unlist(rv.custom$Variablefiltering_query), collapse = " "),
                       nbDeleted = nBefore - nAfter,
                       TotalMainAssay = nrow(rv$dataIn[[length(rv$dataIn)]]))
      
      rv.custom$variable_Filter_SummaryDT <- rbind(rv.custom$variable_Filter_SummaryDT , 
                                                    df)
      
      # Add the parameters values to the new dataset
      par <- rv.custom$Variablefiltering_query
      params(rv$dataIn[[length(rv$dataIn)]]) <- par
      
      #params(rv$dataIn, length(rv$dataIn)) <- reactiveValuesToList(rv.widgets)
      
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Variablefiltering'] <- global$VALIDATED
    })
    
    
    
    
    
    #--------------------------------------------------------
    output$Save <- renderUI({
      wellPanel(
        uiOutput(ns('Save_btn_validate_ui')),
        uiOutput(ns('mod_dl_ui'))
      )
    })
    
    
    output$mod_dl_ui <- renderUI({
      req(config@mode == 'process')
      req(rv$steps.status['Save'] == global$VALIDATED)
        mod_dl_ui(ns('createQuickLink'))
    })
    
    
    output$Save_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Save_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      Magellan::toggleWidget(widget, rv$steps.enabled['Save'])

    })
    
    
    observeEvent(input$Save_btn_validate, {
      
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
      mod_dl_server('createQuickLink', dataIn = reactive({rv$dataIn}))
    })
    
    
    #----------------------------------------------------
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

