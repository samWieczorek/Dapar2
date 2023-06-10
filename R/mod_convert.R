# Module UI


redBtnClass <- "btn-danger"
PrevNextBtnClass <- "btn-info"

optionsBtnClass <- "info"


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

  # This list contains the basic configuration of the process
  config <- Config(
    fullname = 'Convert',
    # Define the type of module
    mode = 'process',
    # List of all steps of the process
    steps = c('Select File'),
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE)
  )

    widgets.default.values <- list(
      Step1_typeOfData = NULL,
      Step1_checkDataLogged = "no",
      Step1_replaceAllZeros = TRUE,
      Step1_software = character(0),
      Step1_XLSsheets = NULL
    )

    ### -------------------------------------------------------------###
    ###                                                             ###
    ### ------------------- MODULE SERVER --------------------------###
    ###                                                             ###
    ### -------------------------------------------------------------###
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Insert necessary code which is hosted by MagellanNTK
        # DO NOT MODIFY THIS LINE
        eval(
          str2expression(
            Get_Workflow_Core_Code(
              w.names = names(widgets.default.values),
              rv.custom.names = names(rv.custom.default.values)
            )
          )
        )
        # >>>
        # >>> START ------------- Code for Description UI---------------
        # >>> 
        
        
        output$Description <- renderUI({
          file <- paste0(config@path_to_md_dir, '/', id, '.md')
          
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
            uiOutput(ns('SelectFile_select3_ui')),
            mod_foo_ui(ns('foo')),
            # Insert validation button
            uiOutput(ns('Step1_btn_validate_ui')),
            
            # Additional code
            plotOutput(ns('showPlot'))
          )
        })
        
        
        # >>> START: Definition of the widgets
        
        
        
        
        rv.custom$mod_foo <- mod_foo_server('foo',
                                            obj = reactive({rv$dataIn}),
                                            reset = reactive({NULL}),
                                            is.enabled = reactive({rv$steps.enabled['Step1']})
        )
        
        
        
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
        
    })
}




#--------------------------------------------------

library(MagellanNTK)

#f <- system.file("module_examples", "extdata/mod_Process1.R", package="MagellanNTK")
#source(file(f), local=TRUE)$value

# ui <- fluidPage(
#   nav_ui('Convert')
# )
# 
# 
# server <- function(input, output){
#   rv <- reactiveValues(
#     dataIn = NULL,
#     dataOut = NULL
#   )
#   
#   observe({
#     rv$dataOut <- nav_server(id = 'Convert')
#   }, priority=1000)
# }
# 
# 
# shinyApp(ui, server)
# 





run_workflow("Convert", dataIn = data.frame())
