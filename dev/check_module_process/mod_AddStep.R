
mod_AddStep_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    uiOutput(ns('form_ui'))
)
}


#' @title xxx
#'
#' @description
#' This module contains the configuration informations for the corresponding pipeline.
#' It is called by the nav_pipeline module of the package Magellan
#'
#' @param id xxx
#'
#' @param dataIn The dataset
#'
#' @param steps.enabled A vector of boolen which has the same length of the steps
#' of the pipeline. xxx
#'
#' @param reset It is a remote command to reset the module. A boolen that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#'
#' @param status xxx
#'
#' @author Samuel Wieczorek
#'
#' @export
#'
#' @import shiny
#'
mod_AddStep_server <- function(id, n, nTotal){

  #' @field global xxxx
  rv <- reactiveValues(
    name = NULL,
    isMandatory = FALSE
  )

  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$form_ui <- renderUI({

      if (n() == nTotal()){
          wellPanel(style="background-color: lightblue",
          fluidRow(
            column(3, shinyjs::disabled(textInput(ns('name'),'Step name', value = 'Save'))),
            column(3, shinyjs::disabled(selectInput(ns('isMandatory'),'Mandatory', choices = c(T), width = '100px'))),
            column(3, selectInput(ns('members'), "Nb widgets", choices = 0:5, width = '100px'))
          ),
          uiOutput(ns('ll_widgets_ui')))
      } else if (n() == 1){
        wellPanel(style="background-color: lightblue",
                  fluidRow(
          column(3, shinyjs::disabled(textInput(ns('name'),'Step name', value = 'Description'))),
          column(3, shinyjs::disabled(selectInput(ns('isMandatory'),'Mandatory', choices = c(T), width = '100px'))),
          column(3, shinyjs::disabled(selectInput(ns('members'), "Nb widgets", choices = 0, width = '100px')))
        ),
        uiOutput(ns('ll_widgets_ui')))
      } else {
        wellPanel(style="background-color: lightgrey",
                  fluidRow(
          column(3, textInput(ns('name'),'Step name')),
          column(3, selectInput(ns('isMandatory'),'Mandatory', choices = c(T, F), width = '100px')),
          column(3, selectInput(ns('members'), "Nb widgets", choices = 0:5, width = '100px'))
        ),
        uiOutput(ns('ll_widgets_ui')))
        }



    })




    output$ll_widgets_ui <- renderUI({
      members <- as.integer(input$members) # default 2
      req(input$members > 0)
      lapply(1:members, function(i) {
        tagList(
          fluidRow(
            column(3, selectInput(ns(paste0('widget_type_', i)),'Type', choices = c('selectInput', 'numericInput', 'radioButtons'), width = '150px')),
            column(3, selectInput(ns(paste0('widget_renderUI_', i)),'In renderUI', choices = c(FALSE, TRUE), width = '100px')),
            column(3, textInput(ns(paste0('widget_name_', i)),'Name')),
            column(3, textInput(ns(paste0('widget_defaultVal_', i)), 'Default value'))
            )
          )
      })
    })

    Get_Value <- function(item){
      #browser()
      members <- as.integer(input$members)
      if(members == 0)
        ll <- NULL
      else
        ll <- lapply(1:members, function(x) {
        req(input[[paste0('widget_', item, '_',x)]])
        input[[paste0('widget_', item, '_',x)]]
        }
        )
      return(ll)
    }


    # Return value of module
    # DO NOT MODIFY THIS PART
    list(name =  reactive({input$name}),
         isMandatory = reactive({as.logical(input$isMandatory)}),
         widgets.type = reactive({Get_Value('type')}),
         widgets.renderUI = reactive({as.logical(Get_Value('renderUI'))}),
         widgets.name = reactive({if (length(Get_Value('name')) > 0) lapply(1:length(Get_Value('name')), function(x){paste0(input$name,'_',  Get_Value('name')[[x]])}) else NULL }),
         widgets.defaultVal = reactive({Get_Value('defaultVal')})
    )

  }
  )
}
