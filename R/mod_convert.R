# Module UI


redBtnClass <- "btn-danger"
PrevNextBtnClass <- "btn-info"
btn_success_color <- "btn-success"
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
mod_Convert_ui <- function(id){
  # ns <- NS(id)
  # tagList()
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
#' @importFrom stats setNames
#'
#' @export
#'
#' @rdname mod_convert
#' 
#' @return NA
#' 
mod_Convert_server <- function(id,
                               dataIn = NULL,
                               steps.enabled = reactive({NULL}),
                               remoteReset = reactive({FALSE})
){

  if (! requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install openxlsx: BiocManager::install('openxlsx')")
  }
  
  global <- list(
    VALIDATED = 1,
    UNDONE = 0,
    SKIPPED = -1
  )

  widgets.default.values <- list(
    selectFile_typeOfData = NULL,
    selectFile_checkDataLogged = 'no',
    selectFile_replaceAllZeros = TRUE,
    selectFile_software = character(0)
  )

  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rv <- reactiveValues(
      dataIn = NULL,
      dataOut = NULL,
      status = NULL,
      reset = NULL,
      steps.enabled = NULL
    )


    dataOut <- reactiveValues(
      trigger = NULL,
      value = NULL
    )

    config <- list(
      name = 'Convert',
      steps = c('Description', 'SelectFile', 'DataID', 'ExpFeatData', 'SamplesMetadata', 'Save'),
      mandatory = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
    )

    rv.widgets <- reactiveValues(
      selectFile_typeOfData = widgets.default.values$selectFile_typeOfData,
      selectFile_checkDataLogged = widgets.default.values$selectFile_checkDataLogged,
      selectFile_replaceAllZeros = widgets.default.values$selectFile_replaceAllZeros,
      selectFile_software = widgets.default.values$selectFile_software,
      selectFile_XLSsheets = widgets.default.values$selectFile_XLSsheets

    )

    #
    # Initialization of the module
    #
    observeEvent(steps.enabled(), ignoreNULL = TRUE, {
      if (is.null(steps.enabled()))
        rv$steps.enabled <- setNames(rep(FALSE, length(config@steps)), config@steps)
      else
        rv$steps.enabled <- steps.enabled()
    })


    observeEvent(remoteReset(), {
      lapply(names(rv.widgets), function(x){
        rv.widgets[[x]] <- widgets.default.values[[x]]
      })
    })


    ###### ------------------- Code for Description (step 0) -------------------------    #####
    output$Description <- renderUI({
      rv$steps.enabled
      
      widget <- actionButton(ns('btn_validate_Description'),
                             paste0('Start ', config@name),
                             class = btn_success_color)
      wellPanel(
        tagList(
          includeMarkdown( system.file('app/md', paste0(config@name, '.md'), 
                                       package = 'DaparToolshed')),
          uiOutput(ns('datasetDescription')),
          Magellan::toggleWidget(widget, rv$steps.enabled['Description'])
        )
      )

      observeEvent(input$btn_validate_Description, ignoreInit = TRUE, ignoreNULL=TRUE, {
        rv$dataIn <- dataIn()
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
      })


    })






    ###-------------------------------------------------------------###
    ### --------------- Code for step 1: selectFile ----------------###
    ###-------------------------------------------------------------###




    output$choose_file_to_import_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      rv.widgets$selectFile_file2convert
      
      widget <- fileInput(ns("file2convert"), "",
                          multiple=FALSE,
                          accept=c(".txt", ".tsv", ".csv",".xls", ".xlsx"))
      
      fluidRow(
        column(width=2, 
               mod_popover_for_help_ui("modulePopover_convertChooseDatafile")),
        column(width = 10,
               Magellan::toggleWidget(widget, rv$steps.enabled['SelectFile'])
        ))
    })


    output$ConvertOptions_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      req(input$file2convert)
      rv.widgets$selectFile_typeOfData
      rv.widgets$selectFile_checkDataLogged
      rv.widgets$selectFile_replaceAllZeros

      widget1 <- radioButtons(ns("typeOfData"),
                              "Is it a peptide or protein dataset ?",
                              choices=c("peptide dataset" = "peptide",
                                        "protein dataset" = "protein"),
                              selected = rv.widgets$selectFile_typeOfData
                              )
      
      widget2 <- radioButtons(ns("checkDataLogged"),
                              "Are your data already log-transformed ?",
                              choices=c("yes (they stay unchanged)" = "yes",
                                        "no (they will be automatically transformed)"="no"),
                              selected = rv.widgets$selectFile_checkDataLogged)
      
      widget3 <- checkboxInput(ns("replaceAllZeros"),
                               "Replace all 0 and NaN by NA",
                               value = rv.widgets$selectFile_replaceAllZeros)
      
      
      tagList(
        Magellan::toggleWidget(widget1, rv$steps.enabled['SelectFile']),
        Magellan::toggleWidget(widget2, rv$steps.enabled['SelectFile']),
        br(),
        Magellan::toggleWidget(widget2, rv$steps.enabled['SelectFile'])
      )

      observeEvent(input$selectFile_typeOfData,{rv.widgets$selectFile_typeOfData <- input$selectFile_typeOfData})
      observeEvent(input$selectFile_checkDataLogged,{rv.widgets$selectFile_checkDataLogged <- input$selectFile_checkDataLogged})
      observeEvent(input$selectFile_replaceAllZeros,{rv.widgets$selectFile_replaceAllZeros <- input$selectFile_replaceAllZeros})

    })


    output$ManageXlsFiles_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      req(input$file2convert)

      if ( GetExtension(input$file2convert$name) %in% c("xls", "xlsx")){
        sheets <- openxlsx::getSheetNames(input$file2convert$datapath)
        widget <- selectInput(ns("selectFile_XLSsheets"), "sheets",
                              choices = as.list(sheets),
                              selected = rv.widgets$selectFile_XLSsheets,
                              width='200px')
        Magellan::toggleWidget(widget2, rv$steps.enabled['SelectFile'])
      }

      observeEvent(input$selectFile_XLSsheets,{
        rv.widgets$selectFile_XLSsheets <- input$selectFile_XLSsheets
        })

    })




    ##
    ## Main renderUI() function for the step SelectFile
    ##
    output$SelectFile <- renderUI({
      name <- 'SelectFile'
      
      widget1 <- div(
        radioButtons(ns("choose_software"), "Software to import from",
                     choices = setNames(nm = c('maxquant', 'proline')),
                     selected = rv.widgets$selectFile_software
        )
      )
      
      widget2 <- div(
        actionButton(ns('btn_validate_SelectFile'),
                     'Perform SelectFile',
                     class = Magellan::btn_success_color)
      )
      
      
      tagList(
        Magellan::toggleWidget(widget1, rv$steps.enabled['SelectFile']),
        uiOutput(ns('choose_file_to_import_ui')),
        uiOutput(ns("ManageXlsFiles_ui")),
        uiOutput(ns("ConvertOptions_ui")),
        Magellan::toggleWidget(widget2, rv$steps.enabled['SelectFile'])
      )
    })


    observeEvent(input$selectFile_software,{
      rv.widgets$selectFile_software <- input$selectFile_software
      })

    observeEvent(input$btn_validate_SelectFile, ignoreInit = TRUE, {
      # Add your stuff code here
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
    })

    ###### ------------------- Code for step 2: DataID -------------------------    #####


    output$DataID <- renderUI({
      name <- 'DataID'
      widget <- actionButton(ns(paste0('btn_validate_', DataID)),
                           'Perform DataID',
                           class = Magellan::btn_success_color)
     Magellan::toggleWidget(widget, rv$steps.enabled['DataID'])
      observeEvent(input$btn_validate_DataID, ignoreInit = TRUE, {
        # Add your stuff code here
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
      })

    })



    ###### ------------------- Code for ExpFeatData -------------------------    #####

    output$ExpFeatData <- renderUI({
      name <- 'ExpFeatData'
      widget <- actionButton(ns(paste0('btn_validate_', ExpFeatData)),
                             'Perform ExpFeatData',
                             class = Magellan::btn_success_color)
      
        Magellan::toggleWidget(widget, rv$steps.enabled['ExpFeatData'])

      observeEvent(input$btn_validate_ExpFeatData, ignoreInit = TRUE, {
        # Add your stuff code here
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
      })

    })



    ###### ------------------- Code for SamplesMetadata -------------------------    #####

    output$SamplesMetadata <- renderUI({
      name <- 'SamplesMetadata'
      widget <- actionButton(ns(paste0('btn_validate_', SamplesMetadata)),
                             'Perform SamplesMetadata',
                             class = btn_success_color)
      
      Magellan::toggleWidget(widget, rv$steps.enabled['SamplesMetadata'])
      

      observeEvent(input$btn_validate_SamplesMetadata, ignoreInit = TRUE, {
        # Add your stuff code here
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
      })

    })



    ###### ------------------- Code for Save -------------------------    #####

    output$Save <- renderUI({
      name <- 'Save'
      widget <- actionButton(ns(paste0('btn_validate_', Save)),
                         'Perform Save',
                         class = btn_success_color)
      Magellan::toggleWidget(widget, rv$steps.enabled['Save'])
      

      observeEvent(input$btn_validate_Save, ignoreInit = TRUE, {
        # Add your stuff code here
        rv$dataIn <- AddItemToDataset(rv$dataIn, config@name)
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
      })


    })

    #------------- Code for validation step ---------------




    # Return value of module
    # DO NOT MODIFY THIS PART
    list(config = reactive({
      config@ll.UI <- setNames(lapply(config@steps,
                                      function(x){
                                        do.call('uiOutput', list(ns(x)))
                                      }),
                               paste0('screen_', config@steps)
      )
      config
    }),
    dataOut = reactive({dataOut})
    #status = reactive({rv$status})
    )


  })
}
