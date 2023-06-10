
redBtnClass <- "btn-danger"
PrevNextBtnClass <- "btn-info"
btn_success_color <- 'info'

optionsBtnClass <- "info"
options(shiny.fullstacktrace = T)

#' @rdname example_module_process1
#' @export
#' 
Convert_conf <- function(){
  # This list contains the basic configuration of the process
  Config(
    fullname = 'Convert',
    # Define the type of module
    mode = 'process',
    # List of all steps of the process
    steps = c('Select File'),
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE)
  )
}



#' @title   mod_choose_pipeline_ui and mod_choose_pipeline_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#'
#' @rdname mod_convert
#'
#' @keywords internal
#' @export
#'
#' @importFrom shiny NS tagList
#' @import sos
#'
#' @return NA
#'
Convert_ui <- function(id) {
  ns <- NS(id)
}




#' Convert Server Function
#'
#' @param id xxx
#' @param dataIn xxx
#' @param steps.enabled xxx
#' @param remoteReset xxx
#'
#' @import QFeatures
#' @importFrom shinyalert shinyalert
#' @importFrom shinyjs disabled
#'
#' @export
#'
#' @rdname mod_convert
#'
#' @return NA
#'
Convert_server <- function(id,
                           dataIn = NULL,
                           steps.enabled = reactive({NULL}),
                           remoteReset = reactive({FALSE}),
                           steps.status = reactive({NULL}),
                           current.pos = reactive({1}),
                           verbose = FALSE
) {
  
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install openxlsx: BiocManager::install('openxlsx')")
  }
  
  
  
  widgets.default.values <- list(
    SelectFile_software = '',
    SelectFile_file = '',
    SelectFile_typeOfData = NULL,
    SelectFile_checkDataLogged = "no",
    SelectFile_replaceAllZeros = TRUE,
    SelectFile_software = character(0),
    SelectFile_XLSsheets = NULL
  )
  
  rv.custom.default.values <- list()
  
  ### -------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ### -------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    core.code <- Get_Workflow_Core_Code(
      mode = 'process',
      name = id,
      w.names = names(widgets.default.values),
      rv.custom.names = names(rv.custom.default.values)
    )
    
    eval(str2expression(core.code))
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      md.file <- paste0(id, '.md')
      path <- system.file('extdata/workflow/PipelineA/md', package='MagellanNTK')
      file <- file.path(path, md.file)
      
      tagList(
        # In this example, the md file is found in the extdata/module_examples directory
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
      # Insert your own code to visualize some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"),
                             "Start",
                             class = btn_success_color)
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- global$VALIDATED
    })
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$SelectFile <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        uiOutput(ns('SelectFile_software_ui')),
        uiOutput(ns('SelectFile_file_ui')),
        uiOutput(ns('SelectFile_ManageXlsFiles_ui')),
        uiOutput(ns('SelectFile_typeOfData_ui')),
        uiOutput(ns('SelectFile_checkDataLogged_ui')),
        uiOutput(ns('SelectFile_replaceAllZeros_ui')),
        
        # Insert validation button
        uiOutput(ns('SelectFile_btn_validate_ui')),
        
      )
    })
    
    
    # >>> START: Definition of the widgets
    
    
    
    output$SelectFile_software_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_software"), 
                             "Software to import from",
                             choices = setNames(nm = DAPAR::GetSoftAvailables()),
                             selected = character(0)
      )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    # This part must be customized by the developer of a new module
    output$SelectFile_file_ui <- renderUI({
      req(input$choose_software)
      fluidRow(
        column(width = 2,
               popover_for_help_ui("modulePopover_convertChooseDatafile")
        ),
        column(width = 10,
               widget <- fileInput(ns("file"), "",
                                   multiple = FALSE,
                                   accept = c(".txt", ".tsv", ".csv", ".xls", ".xlsx")
               )
        )
      )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    fileExt.ok <- reactive({
      req(input$file1$name)
      authorizedExts <- c("txt", "csv", "tsv", "xls", "xlsx")
      ext <- GetExtension(input$file1$name)
      !is.na(match(ext, authorizedExts))
    })
    
    
    
    output$SelectFile_ManageXlsFiles_ui <- renderUI({
      req(input$choose_software)
      req(input$file1)
      
      req(GetExtension(input$file1$name) %in% c("xls", "xlsx"))
      
      tryCatch({   
        sheets <- listSheets(input$file1$datapath)
        widget <- selectInput(ns("XLSsheets"), 
                              "sheets", 
                              choices = as.list(sheets), 
                              width = "200px")
      },
      warning = function(w) {
        shinyjs::info(conditionMessage(w))
        return(NULL)
      },
      error = function(e) {
        shinyjs::info(conditionMessage(e))
        return(NULL)
      },
      finally = {
        # cleanup-code
      }
      )
      toggleWidget(widget, rv$steps.enabled['SelectFile'])
    })
    
    
    
    output$SelectFile_typeOfData_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_typeOfData"), 
                             "Is it a peptide or protein dataset ?",
                             choices = c(
                               "peptide dataset" = "peptide",
                               "protein dataset" = "protein"
                             )
                             )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    output$SelectFile_checkDataLogged_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_checkDataLogged"), 
                             "Are your data already log-transformed ?",
                             choices = c(
                               "yes (they stay unchanged)" = "yes",
                               "no (they wil be automatically transformed)" = "no"
                             ),
                             selected = "no"
                             )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    output$SelectFile_replaceAllZeros_ui <- renderUI({
      widget <- checkboxInput(ns("SelectFile_replaceAllZeros"), 
                             "Replacve all 0 and NaN by NA",
                             value = TRUE
                             )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    
    output$SelectFile_btn_validate_ui <- renderUI({
      widget <-  actionButton(ns("SelectFile_btn_validate"),
                              "Perform",
                              class = btn_success_color)
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
      
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$SelectFile_btn_validate, {
      # Do some stuff
      rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
                                          dataset = rnorm(1:5),
                                          name = paste0('temp_',id))
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['SelectFile'] <- global$VALIDATED
    })
    
    
    # <<< END ------------- Code for step 1 UI---------------
    
    # >>> START ------------- Code for step 3 UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
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
      toggleWidget(actionButton(ns("Save_btn_validate"), "Save",
                                class = btn_success_color),
                   rv$steps.enabled['Save']
      )
    })
    observeEvent(input$Save_btn_validate, {
      # Do some stuff
      rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
                                          dataset = rnorm(1:5),
                                          name = id)
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
      mod_dl_server('createQuickLink', dataIn = reactive({rv$dataIn}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  })
}




#--------------------------------------------------

library(MagellanNTK)

run_workflow("Convert", dataIn = data.frame())
