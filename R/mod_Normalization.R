#' @title Shiny example process module.
#'
#' @description
#' This module contains the configuration informations for the corresponding pipeline.
#' It is called by the nav_pipeline module of the package Magellan
#' 
#' The name of the server and ui functions are formatted with keywords separated by '_', as follows:
#' * first string `mod`: indicates that it is a Shiny module
#' * `pipeline name` is the name of the pipeline to which the process belongs
#' * `process name` is the name of the process itself
#' 
#' This convention is important because Magellan call the different
#' server and ui functions by building dynamically their name.
#' 
#' In this example, `mod_Normalization_ui()` and `mod_Normalization_server()` define
#' the code for the process `Normalization`.
#' 
#' @param id xxx
#' 
#' @rdname example_module_Normalization
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#'
mod_Normalization_ui <- function(id){
  ns <- NS(id)
}


#' @param id xxx
#'
#' @param dataIn The dataset
#'
#' @param steps.enabled A vector of boolean which has the same length of the steps
#' of the pipeline. This information is used to enable/disable the widgets. It is not
#' a communication variable between the caller and this module, thus there is no
#' corresponding output variable
#'
#' @param remoteReset It is a remote command to reset the module. A boolean that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#' 
#' @param steps.status xxx
#' 
#' @param current.pos xxx
#'
#' @rdname example_module_Normalization
#' 
#' @importFrom stats setNames rnorm
#' 
#' @export
#' 
mod_Normalization_server <- function(id,
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
    steps = c('Description', 'Normalization', 'Save'),
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, TRUE, TRUE),
    
    path_to_md_dir = system.file('module_examples/md/', package='Magellan')
  )
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    Normalization_method = "None",
    Normalization_type = "overall",
    Normalization_varReduction = FALSE,
    Normalization_quantile = 0.15,
    Normalization_spanLOESS = 0.7,
    Normalization_trackFromBoxplot = NULL,
    Normalization_selectProt = NULL, 
    Normalization_resetTracking = FALSE,
    Normalization_sync = FALSE
  )
  
  
  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    eval(str2expression(Get_Worflow_Core_Code(
      w.names = names(widgets.default.values)
    )))
    
    rv.custom <- reactiveValues()
    
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
      toggleWidget(rv$steps.enabled['Description'], widget)
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- global$VALIDATED
    })
    
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$Normalization <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by Magellan
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        uiOutput(ns('Step1_btn1_ui')),
        uiOutput(ns('Step1_select1_ui')),
        uiOutput(ns('Step1_select2_ui')),
        uiOutput(ns('Step1_select3_ui')),
        # Insert validation button
        uiOutput(ns('Step1_btn_validate_ui')),
        
        # Additional code
        plotOutput(ns('showPlot'))
        
        
        
        
      )
      tagList(
        div(
          div(
            style="display:inline-block; vertical-align: middle; padding-right: 20px;",
            uiOutput(ns('Normalization_method_ui'))
          ),
          div(
            style="display:inline-block; vertical-align: middle; padding-right: 20px;",
            hidden(selectInput("normalization.type", "Normalization type",  
                               choices = setNames( nm = c("overall", "within conditions")), 
                               selected = rv$widgets$normalization$type,
                               width='150px'))
          ),
          div(
            style="display:inline-block; vertical-align: middle; padding-right: 20px;",
            hidden(textInput("spanLOESS", "Span",value = rv$widgets$normalization$spanLOESS, width='100px')),
            module_Not_a_numericUI("test_spanLOESS"),
            uiOutput("choose_normalizationQuantile"),
            uiOutput("choose_normalizationScaling")
          ),
          hidden(
            div(id = 'DivMasterProtSelection',
                style="display:inline-block; vertical-align: middle; padding-right: 20px;",
                mod_plots_tracking_ui('master_tracking'),
                checkboxInput("SyncForNorm", "Synchronise with selection above", value=FALSE)
            )
          ),
          
          div(
            style="display:inline-block; vertical-align: middle; padding-right: 20px;",
            hidden(actionButton("perform.normalization", "Perform normalization", class = actionBtnClass, width="170px"))
          )
        ),
        uiOutput("helpForNormalizationMethods"),
        tags$hr(),
        fluidRow(
          column(width=4, moduleDensityplotUI("densityPlot_Norm")),
          column(width=4,
                 withProgress(message = 'Building plot',detail = '', value = 0, {
                   mod_plots_intensity_ui("boxPlot_Norm")
                 })),
          column(width=4,withProgress(message = 'Building plot',detail = '', value = 0, {
            highchartOutput("viewComparisonNorm_HC")
          })
          )
        )
      )
    })
    
    
    })
    
    
    # >>> START: Definition of the widgets
    
    output$Step1_btn1_ui <- renderUI({
      
      
    })
    
    output$Normalization_method_ui <- renderUI({
      widget <- selectInput("Normalization_method","Normalization method", 
                            choices = normMethods, 
                            selected = rv.widgets$Normalization_method,
                            width='200px')
      toggleWidget(rv$steps.enabled['Normalization'], widget )
    })
    
    # This part must be customized by the developer of a new module
    output$Step1_select1_ui <- renderUI({
      widget <- selectInput(ns('Step1_select1'),
                            'Select 1 in renderUI',
                            choices = 1:4,
                            selected = rv.widgets$Step1_select1,
                            width = '150px')
      toggleWidget(rv$steps.enabled['Step1'], widget )
    })
    
    
    output$Step1_select2_ui <- renderUI({
      widget <- selectInput(ns('Step1_select2'),
                            'Select 2 in renderUI',
                            choices = 1:4,
                            selected = rv.widgets$Step1_select2,
                            width = '150px')
      toggleWidget(rv$steps.enabled['Step1'], widget )
    })
    
    
    output$Step1_select3_ui <- renderUI({
      widget <- selectInput(ns('Step1_select3'),
                            'Select 1 in renderUI',
                            choices = 1:4,
                            selected = rv.widgets$Step1_select3,
                            width = '150px')
      toggleWidget(rv$steps.enabled['Step1'], widget )
    })
    
    
    
    output$Step1_btn_validate_ui <- renderUI({
      widget <-  actionButton(ns("Step1_btn_validate"),
                              "Perform",
                              class = btn_success_color)
      toggleWidget(rv$steps.enabled['Step1'], widget )
      
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$Step1_btn_validate, {
      # Do some stuff
      
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Step1'] <- global$VALIDATED
    })
    
    
    output$showPlot <- renderPlot({
      plot(as.matrix(dataIn()[[1]]))
    })
    # <<< END ------------- Code for step 1 UI---------------
    
    
    # >>> START ------------- Code for step 2 UI---------------
    
    output$Step2 <- renderUI({
      wellPanel(
        # Two examples of widgets in a renderUI() function
        uiOutput(ns('Step2_select1_ui')),
        uiOutput(ns('Step2_select2_ui')),
        
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Step2_btn_validate_ui'))
      )
    })
    
    
    output$Step2_select1_ui <- renderUI({
      widget <- selectInput(ns('Step2_select1'),
                            'Select 1 in renderUI',
                            choices = 1:4,
                            selected = rv.widgets$Step2_select1,
                            width = '150px')
      toggleWidget(rv$steps.enabled['Step2'], widget )
    })
    
    output$Step2_select2_ui <- renderUI({
      widget <- selectInput(ns('Step2_select2'),
                            'Select 1 in renderUI',
                            choices = 1:4,
                            selected = rv.widgets$Step2_select2,
                            width = '150px')
      toggleWidget(rv$steps.enabled['Step2'], widget )
    })
    
    output$Step2_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Step2_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      toggleWidget(rv$steps.enabled['Step2'], widget )
    })
    
    observeEvent(input$Step2_btn_validate, {
      # Do some stuff
      
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Step2'] <- global$VALIDATED
    })
    
    # <<< END ------------- Code for step 2 UI---------------
    
    
    # >>> START ------------- Code for step 3 UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Save_btn_validate_ui'))
      )
    })
    
    output$Save_btn_validate_ui <- renderUI({
      tagList(
        toggleWidget(rv$steps.enabled['Save'], 
                     actionButton(ns("Save_btn_validate"), "Save",
                                  class = btn_success_color)
        ),
        if (config$mode == 'process' && rv$steps.status['Save'] == global$VALIDATED) {
          mod_Save_Dataset_ui(ns('createQuickLink'))
        }
      )
      
    })
    observeEvent(input$Save_btn_validate, {
      # Do some stuff
      rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
                                          dataset = rnorm(1:5),
                                          name = id)
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
      mod_Save_Dataset_server('createQuickLink', dataIn = reactive({rv$dataIn}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    
    
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}
