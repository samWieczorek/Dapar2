
#' proc1 Server Function
#'
#' @noRd
#'
#' @import QFeatures
#' @importFrom shinyalert shinyalert
#'
mod_proc1_server <- function(id,
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
# Define default selected values for widgets
widgets.default.values <- list(
Step1_a = a,
Step1_b = b
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
    name = 'proc1',
    steps = c('Description','Step1','Save'),
    mandatory = c(TRUE,TRUE,TRUE)
  )
rv.widgets <- reactiveValues(
Step1_a =  widgets.default.values$Step1_a,
Step1_b =  widgets.default.values$Step1_b
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


###### ------------------- Code for Step1 -------------------------    #####

observeEvent(input$Step1_a,{rv.widgets$Step1_a <- input$Step1_a})

observeEvent(input$Step1_b,{rv.widgets$Step1_b <- input$Step1_b})

output$Step1_a_UI <- renderUI({
  rv.widgets$Step1_a
  if (rv$steps.enabled['Step1'])
    selectInput(ns('Step1_a'), 'Step1_a in renderUI',
                choices = 1:4,
                selected = rv.widgets$Step1_a,
                width = '150px')
  else
    shinyjs::disabled(
      selectInput(ns('Step1_a'), 'Step1_a in renderUI',
                  choices = 1:4,
                  selected = rv.widgets$Step1_a,
                  width = '150px')
      )
})

output$Step1 <- renderUI({
  name <- 'Step1'
  wellPanel(id = ns('foo'),
    actionButton(ns('exampleBtn_Step1'), 'Example button for Step1'),
    tagList(

   uiOutput(ns('Step1_a_UI'))
   ,
div(
if (rv$steps.enabled['Step1'])
                selectInput(ns('Step1_b'), 'Select Step1_b',
                          choices = 1:3,
                          selected = rv.widgets$Step1_b,
                          width = '150px')
              else
                shinyjs::disabled(
                  selectInput(ns('Step1_b'), 'Select Step1_b',
                              choices = 1:5,
                              selected = rv.widgets$Step1_b,
                              width = '150px')
                )
)
,
div(
if (rv$steps.enabled['Step1'])
                actionButton(ns(paste0('btn_validate_', Step1)),
                             'Perform Step1',
                             class = btn_success_color)
              else
                shinyjs::disabled(
                  actionButton(ns(paste0('btn_validate_', Step1)),
                               'Perform Step1',
                               class = btn_success_color)
                  )
)

    )
  )
})

observeEvent(input$btn_validate_Step1, ignoreInit = T, {
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

