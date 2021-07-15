
mod_Protein_Description_ui <- function(id){
  ns <- NS(id)
}


header <- "btn_style <- 'display:inline-block; vertical-align: middle; padding: 7px'"

ui <- "
#' #name# UI Function
#'
#' @description A shiny Module.
#'
#' @param id, input, output, session  #Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import shinyjs
#' @importFrom shinyalert useShinyalert
#'
mod_#name#_ui <- function(id){
  ns <- NS(id)
}
"



start_server <- "
#' #name# Server Function
#'
#' @noRd
#'
#' @import QFeatures
#' @importFrom shinyalert shinyalert
#'
mod_#name#_server <- function(id,
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
"


moduleServerFunc_template <- "
 ###-------------------------------------------------------------###
 ###                                                             ###
 ### ------------------- MODULE SERVER --------------------------###
 ###                                                             ###
 ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    # Do not modify this part
    ns <- session$ns

    # Do not modify this part
    rv <- reactiveValues(
      dataIn = NULL,
      dataOut = NULL,
      status = NULL,
      reset = NULL,
      steps.enabled = NULL
    )

    # Do not modify this part
    #' @field dataOut xxx
    dataOut <- reactiveValues(
      trigger = NULL,
      value = NULL
    )
"

config_template <- "
config <- list(
    name = '#name#',
    steps = #steps#,
    mandatory = #mandatory#
  )
"


chunk_init_module_template <- "
 ###-------------------------------------------------------------###
 ### ------------ Initialization of the module ------------------###
 ###-------------------------------------------------------------###
  "

chunk_observeEvent_steps_enabled <-"
    observeEvent(steps.enabled(), ignoreNULL = TRUE, {
      if (is.null(steps.enabled()))
        rv$steps.enabled <- setNames(rep(FALSE, rv.process$length), rv.process$config$steps)
      else
        rv$steps.enabled <- steps.enabled()
    })
"

chunk_observeEvent_remoteReset <-"
    observeEvent(remoteReset(), {
      lapply(names(rv.widgets), function(x){
        rv.widgets[[x]] <- widgets.default.values[[x]]
      })
    })
"


output_Description_template <- "
###-------------------------------------------------------------###
### -------------- Code for step #n#: Description --------------###
###-------------------------------------------------------------###
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

"
write_comment_for_code <- "
###-------------------------------------------------------------###
### -------------- Code for step #n#: #step.name# --------------###
###-------------------------------------------------------------###
"
observer_for_widget_comment <- "
##
## Code for listeners of reactive varaibles in the step
##
"
observer_for_widget_template <- "
observeEvent(input$#name#,{rv.widgets$#name# <- input$#name#})
"


widget_renderUI_template <-"
##
## Widget embedded in an renderUI function
##

output$#name#_UI <- renderUI({
  rv.widgets$#name#
  div(style='display:inline-block; vertical-align: middle; padding-right: 40px;',
       if (rv$steps.enabled['#step#'])
    selectInput(ns('#name#'), '#name# in renderUI',
                choices = 1:4,
                selected = rv.widgets$#name#,
                width = '150px')
  else
    shinyjs::disabled(
      selectInput(ns('#name#'), '#name# in renderUI',
                  choices = 1:4,
                  selected = rv.widgets$#name#,
                  width = '150px')
      )
  )
})
"


code_for_global_step_renderUI_header <- "
##
## Main renderUI() function for the step #step.name#
##
output$#step.name# <- renderUI({
  name <- '#step.name#'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_#step.name#'), 'Example button for #step.name#'),
    tagList(
"

code_for_global_step_renderUI_call_widget_renderUI <- "
   div(style='display:inline-block; vertical-align: middle; padding-right: 40px;',
       uiOutput(ns('#widget.name#_UI'))
   )
   "

code_for_global_step_renderUI_call_widget <- "
##
## Widget included in the global renderUI() function of the step
##
div(
if (rv$steps.enabled['#step.name#'])
                selectInput(ns('#widget.name#'), 'Select #widget.name#',
                          choices = 1:3,
                          selected = rv.widgets$#widget.name#,
                          width = '150px')
              else
                shinyjs::disabled(
                  selectInput(ns('#widget.name#'), 'Select #widget.name#',
                              choices = 1:5,
                              selected = rv.widgets$#widget.name#,
                              width = '150px')
                )
)
"

code_for_global_step_renderUI_call_validate_btn <- "
div(
if (rv$steps.enabled['#step.name#'])
                actionButton(ns(paste0('btn_validate_', #step.name#)),
                             'Perform #step.name#',
                             class = btn_success_color)
              else
                shinyjs::disabled(
                  actionButton(ns(paste0('btn_validate_', #step.name#)),
                               'Perform #step.name#',
                               class = btn_success_color)
                  )
)
"

code_for_global_step_renderUI_end <- "
    )
  )
})
"

chunk_for_observe_validation_btn_generic_step <- "
observeEvent(input$btn_validate_#step.name#, ignoreInit = T, {
  # Add your stuff code here

  # Do not modify this part
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})
"

chunk_for_observe_validation_btn_last_step <- "
#------------- Code for validation step ---------------

observeEvent(input$btn_validate_#step.name#, ignoreInit = T, {
  # Add your stuff code here
  rv$dataIn <- AddItemToDataset(rv$dataIn, config$name)
  dataOut$trigger <- Send_Result_to_Caller(rv$dataIn)$trigger
  dataOut$value <- Send_Result_to_Caller(rv$dataIn)$value
})
"

chunk_for_return_server <-  "
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
"


chunk_for_end_module <- "
})
}

"
