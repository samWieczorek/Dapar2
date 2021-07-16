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
#'
#' @export
#'
#' @rdname mod_convert
#' 
mod_Convert_server <- function(id,
                               dataIn = NULL,
                               steps.enabled = reactive({NULL}),
                               remoteReset = reactive({FALSE})
){

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
      steps = c('Description','SelectFile','DataID','ExpFeatData','SamplesMetadata','Save'),
      mandatory = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
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
        rv$steps.enabled <- setNames(rep(FALSE, rv.process$length), rv.process$config$steps)
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
      #browser()
      wellPanel(
        tagList(
          includeMarkdown( system.file('app/md', paste0(config$name, '.md'), package='Prostar')),
          uiOutput(ns('datasetDescription')),
          if (isTRUE(rv$steps.enabled['Description'])  )
            actionButton(ns('btn_validate_Description'),
                         paste0('Start ', config$name),
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns('btn_validate_Description'),
                           paste0('Start ', config$name),
                           class = btn_success_color)
            )
        )
      )

      # observeEvent(input$btn_validate_Description, ignoreInit = T, ignoreNULL=T, {
      #   rv$dataIn <- dataIn()
      #   dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
      #   dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
      # })


    })






    ###-------------------------------------------------------------###
    ### --------------- Code for step 1: selectFile ----------------###
    ###-------------------------------------------------------------###




    output$choose_file_to_import_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      rv.widgets$selectFile_file2convert
      fluidRow(
        column(width=2, modulePopoverUI("modulePopover_convertChooseDatafile")),
        column(width = 10, if (rv$steps.enabled['selectFile'])
          fileInput(ns("file2convert"), "",
                    multiple=FALSE,
                    accept=c(".txt", ".tsv", ".csv",".xls", ".xlsx"))
          else
            shinyjs::disabled(fileInput(ns("file2convert"), "",
                                        multiple=FALSE,
                                        accept=c(".txt", ".tsv", ".csv",".xls", ".xlsx")))
        ))
    })


    output$ConvertOptions_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      req(input$file2convert)
      rv.widgets$selectFile_typeOfData
      rv.widgets$selectFile_checkDataLogged
      rv.widgets$selectFile_replaceAllZeros

      tagList(
        if (rv$steps.enabled['selectFile'])
          radioButtons(ns("typeOfData"),
                       "Is it a peptide or protein dataset ?",
                       choices=c("peptide dataset" = "peptide",
                                 "protein dataset" = "protein"),
                       selected = rv.widgets$selectFile_typeOfData
          )
        else
          shinyjs::disabled(radioButtons(ns("typeOfData"),
                                         "Is it a peptide or protein dataset ?",
                                         choices=c("peptide dataset" = "peptide",
                                                   "protein dataset" = "protein"),
                                         selected = rv.widgets$selectFile_typeOfData
          )
          )

        , if (rv$steps.enabled['selectFile'])
          radioButtons(ns("checkDataLogged"),
                       "Are your data already log-transformed ?",
                       choices=c("yes (they stay unchanged)" = "yes",
                                 "no (they will be automatically transformed)"="no"),
                       selected = rv.widgets$selectFile_checkDataLogged)
        else
          shinyjs::disabled(radioButtons(ns("checkDataLogged"),
                                         "Are your data already log-transformed ?",
                                         choices=c("yes (they stay unchanged)" = "yes",
                                                   "no (they will be automatically transformed)"="no"),
                                         selected = rv.widgets$selectFile_checkDataLogged)
          )
        ,br()
        ,if (rv$steps.enabled['selectFile'])
          checkboxInput(ns("replaceAllZeros"),
                        "Replace all 0 and NaN by NA",
                        value = rv.widgets$selectFile_replaceAllZeros)
        else
          shinyjs::disabled(checkboxInput(ns("replaceAllZeros"),
                                          "Replace all 0 and NaN by NA",
                                          value = rv.widgets$selectFile_replaceAllZeros)
          )
      )

      observeEvent(input$selectFile_typeOfData,{rv.widgets$selectFile_typeOfData <- input$selectFile_typeOfData})
      observeEvent(input$selectFile_checkDataLogged,{rv.widgets$selectFile_checkDataLogged <- input$selectFile_checkDataLogged})
      observeEvent(input$selectFile_replaceAllZeros,{rv.widgets$selectFile_replaceAllZeros <- input$selectFile_replaceAllZeros})

    })


    output$ManageXlsFiles_ui <- renderUI({
      req(rv.widgets$selectFile_software)
      req(input$file2convert)

      if ( GetExtension(input$file2convert$name) %in% c("xls", "xlsx")){
        sheets <- listSheets(input$file2convert$datapath)
        if (rv$steps.enabled['selectFile'])
          selectInput(ns("selectFile_XLSsheets"), "sheets",
                      choices = as.list(sheets),
                      selected = rv.widgets$selectFile_XLSsheets,
                      width='200px')
        else
          shinyjs::disabled(
            selectInput(ns("selectFile_XLSsheets"), "sheets",
                        choices = as.list(sheets),
                        selected = rv.widgets$selectFile_XLSsheets,
                        width='200px')
          )
      }

      observeEvent(input$selectFile_XLSsheets,{rv.widgets$selectFile_XLSsheets <- input$selectFile_XLSsheets})

    })




    ##
    ## Main renderUI() function for the step SelectFile
    ##
    output$SelectFile <- renderUI({
      name <- 'SelectFile'
      tagList(
        div(style="display:inline-block; vertical-align: middle; padding-right: 40px;",
            if (rv$steps.enabled['SelectFile'])
              radioButtons(ns("choose_software"), "Software to import from",
                           choices = setNames(nm = c('maxquant', 'proline')),
                           selected = rv.widgets$selectFile_software
              )
            else
              shinyjs::disabled(radioButtons(ns("choose_software"), "Software to import from",
                                             choices = setNames(nm = c('maxquant', 'proline')),
                                             selected =  rv.widgets$selectFile_software
              )
              )
        ),
        uiOutput(ns('choose_file_to_import_ui')),
        uiOutput(ns("ManageXlsFiles_ui")),
        uiOutput(ns("ConvertOptions_ui")),
        div(
          if (rv$steps.enabled['SelectFile'])
            actionButton(ns('btn_validate_SelectFile'),
                         'Perform SelectFile',
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns('btn_validate_SelectFile'),
                           'Perform SelectFile',
                           class = btn_success_color)
            )
        )
      )
    })


    observeEvent(input$selectFile_software,{rv.widgets$selectFile_software <- input$selectFile_software})

    observeEvent(input$btn_validate_SelectFile, ignoreInit = T, {
      # Add your stuff code here
      dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
      dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
    })

    ###### ------------------- Code for step 2: DataID -------------------------    #####


    output$DataID <- renderUI({
      name <- 'DataID'
      tagList(

        div(
          if (rv$steps.enabled['DataID'])
            actionButton(ns(paste0('btn_validate_', DataID)),
                         'Perform DataID',
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns(paste0('btn_validate_', DataID)),
                           'Perform DataID',
                           class = btn_success_color)
            )
        )
      )

      observeEvent(input$btn_validate_DataID, ignoreInit = T, {
        # Add your stuff code here
        dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
        dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
      })

    })



    ###### ------------------- Code for ExpFeatData -------------------------    #####

    output$ExpFeatData <- renderUI({
      name <- 'ExpFeatData'
      tagList(

        div(
          if (rv$steps.enabled['ExpFeatData'])
            actionButton(ns(paste0('btn_validate_', ExpFeatData)),
                         'Perform ExpFeatData',
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns(paste0('btn_validate_', ExpFeatData)),
                           'Perform ExpFeatData',
                           class = btn_success_color)
            )
        )
      )

      observeEvent(input$btn_validate_ExpFeatData, ignoreInit = T, {
        # Add your stuff code here
        dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
        dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
      })

    })



    ###### ------------------- Code for SamplesMetadata -------------------------    #####

    output$SamplesMetadata <- renderUI({
      name <- 'SamplesMetadata'
      tagList(

        div(
          if (rv$steps.enabled['SamplesMetadata'])
            actionButton(ns(paste0('btn_validate_', SamplesMetadata)),
                         'Perform SamplesMetadata',
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns(paste0('btn_validate_', SamplesMetadata)),
                           'Perform SamplesMetadata',
                           class = btn_success_color)
            )
        )
      )

      observeEvent(input$btn_validate_SamplesMetadata, ignoreInit = T, {
        # Add your stuff code here
        dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
        dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
      })

    })



    ###### ------------------- Code for Save -------------------------    #####

    output$Save <- renderUI({
      name <- 'Save'
      tagList(
        div(
          if (rv$steps.enabled['Save'])
            actionButton(ns(paste0('btn_validate_', Save)),
                         'Perform Save',
                         class = btn_success_color)
          else
            shinyjs::disabled(
              actionButton(ns(paste0('btn_validate_', Save)),
                           'Perform Save',
                           class = btn_success_color)
            )
        )
      )

      observeEvent(input$btn_validate_Save, ignoreInit = T, {
        # Add your stuff code here
        rv$dataIn <- AddItemToDataset(rv$dataIn, config$name)
        dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
        dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
      })


    })

    #------------- Code for validation step ---------------




    # Return value of module
    # DO NOT MODIFY THIS PART
    list(config = reactive({
      config$ll.UI <- setNames(lapply(config$steps,
                                      function(x){
                                        do.call('uiOutput', list(ns(x)))
                                      }),
                               paste0('screen_', config$steps)
      )
      config
    }),
    dataOut = reactive({dataOut})
    #status = reactive({rv$status})
    )


  })
}
