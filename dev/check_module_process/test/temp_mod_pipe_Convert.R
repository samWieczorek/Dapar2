

#' @export
#'
mod_Convert_ui <- function(id){
  ns <- NS(id)
}




#' Convert Server Function
#'
#' @noRd
#'
#' @import QFeatures
#' @importFrom shinyalert shinyalert
#'
mod_Convert_server <- function(id,
                            dataIn = NULL,
                            steps.enabled = reactive({NULL}),
                            remoteReset = reactive({FALSE})
                            ){

  #' @field global xxxx
  global <- list(
    VALIDATED = 1,
    UNDONE = 0,
    SKIPPED = -1
  )



 ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rv.widgets <- reactiveValues()

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

  # Initialization of the module

###### ------------------- Code for Description (step 0) -------------------------    #####
output$Description <- renderUI({
  rv$steps.enabled
  #browser()
  wellPanel(
    tagList(
      includeMarkdown( system.file('app/md', paste0(config$name, '.md'), package='Magellan')),
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
})



observeEvent(input$btn_validate_Description, ignoreInit = T, ignoreNULL=T, {
  rv$dataIn <- dataIn()
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})


###### ------------------- Code for SelectFile -------------------------    #####

output$SelectFile <- renderUI({
  name <- 'SelectFile'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_SelectFile'), 'Example button for SelectFile'),
    tagList(

div(
if (rv$steps.enabled['SelectFile'])
                actionButton(ns(paste0('btn_validate_', SelectFile)),
                             'Perform SelectFile',
                             class = btn_success_color)
              else
                shinyjs::disabled(
                  actionButton(ns(paste0('btn_validate_', SelectFile)),
                               'Perform SelectFile',
                               class = btn_success_color)
                  )
)

    )
  )
})

observeEvent(input$btn_validate_SelectFile, ignoreInit = T, {
  # Add your stuff code here
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})

###### ------------------- Code for DataID -------------------------    #####

output$DataID <- renderUI({
  name <- 'DataID'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_DataID'), 'Example button for DataID'),
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
  )
})

observeEvent(input$btn_validate_DataID, ignoreInit = T, {
  # Add your stuff code here
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})

###### ------------------- Code for ExpFeatData -------------------------    #####

output$ExpFeatData <- renderUI({
  name <- 'ExpFeatData'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_ExpFeatData'), 'Example button for ExpFeatData'),
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
  )
})

observeEvent(input$btn_validate_ExpFeatData, ignoreInit = T, {
  # Add your stuff code here
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})

###### ------------------- Code for SamplesMetadata -------------------------    #####

output$SamplesMetadata <- renderUI({
  name <- 'SamplesMetadata'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_SamplesMetadata'), 'Example button for SamplesMetadata'),
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
  )
})

observeEvent(input$btn_validate_SamplesMetadata, ignoreInit = T, {
  # Add your stuff code here
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})

###### ------------------- Code for Save -------------------------    #####

output$Save <- renderUI({
  name <- 'Save'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_Save'), 'Example button for Save'),
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
  )
})

#------------- Code for validation step ---------------

observeEvent(input$btn_validate_Save, ignoreInit = T, {
  # Add your stuff code here
  rv$dataIn <- AddItemToDataset(rv$dataIn, config$name)
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})

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

