library(highcharter)
library(shinyjs)
library(DT)
library(tibble)
library(QFeatures)
library(Magellan)
library(DaparToolshed)


source(file.path('../../../R', 'global.R'), local=TRUE)$value
source(file.path('../../../R', 'mod_Protein_Description.R'), local=TRUE)$value
#source(file.path('../../../R', 'mod_plots_tracking.R'), local=TRUE)$value
source(file.path('../../../R', 'mod_settings.R'), local=TRUE)$value
#source(file.path('../../../R', 'mod_plots_density.R'), local=TRUE)$value
#source(file.path('../../../R', 'mod_plots_intensity.R'), local=TRUE)$value
source(file.path("../../../R", "mod_popover_for_help.R"), local = TRUE)$value
#source(file.path('../../../R', 'mod_infos_dataset.R'), local=TRUE)$value
#source(file.path('../../../R', 'mod_observe_dynamic_colourPicker_input.R'), local=TRUE)$value




options(shiny.fullstacktrace = TRUE)



#--------------------------------------------



mod_test_process_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('UI')),
    wellPanel(title = 'foo',
              tagList(
                h3('Valler'),
                uiOutput(ns('show_Debug_Infos'))
              )
    )
  )
}


mod_test_process_server <- function(id){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    utils::data(Exp1_R25_prot, package='DAPARdata2')
    
    obj <- NULL
    obj <- Exp1_R25_prot
    
    rv <- reactiveValues(
      dataIn = Exp1_R25_prot,
      dataOut = NULL
    )

    observe({
       
      rv$dataOut <- mod_nav_process_server(id = 'Protein_Description',
                                           dataIn = reactive({rv$dataIn})
                                           )
      
      
      
      observeEvent(rv$dataOut$dataOut()$trigger, {
        print('totototo')
        print(names(rv$dataOut$dataOut()$value))
      })
    }, priority=1000)
    
    
    output$UI <- renderUI({
      mod_nav_process_ui(ns('Protein_Description'))
    })
    
    
    
    #--------------------------------------------------------------------
    
    output$show_Debug_Infos <- renderUI({
      fluidRow(
        column(width=2,
               tags$b(h4(style = 'color: blue;', "Data In")),
               uiOutput('show_rv_dataIn')),
        column(width=2,
               tags$b(h4(style = 'color: blue;', "Data Out")),
               uiOutput('show_rv_dataOut'))
      )
    })
    
    ###########---------------------------#################
    output$show_rv_dataIn <- renderUI({
      req(rv$dataIn)
      tagList(
        lapply(names(rv$dataIn), function(x){tags$p(x)})
      )
    })
    
    output$show_rv_dataOut <- renderUI({
      req(rv$dataOut)
      tagList(
        lapply(names(rv$dataOut$value), function(x){tags$p(x)})
      )
    })
    
  })
}


#----------------------------------------------------------------------
ui <- fluidPage(
  mod_test_process_ui('test_mod_process')
)


#----------------------------------------------------------------------
server <- function(input, output){
  mod_test_process_server('test_mod_process')
}
shinyApp(ui, server)


# shinyApp(ui = mod_test_process_ui('test_mod_process'),
#          server = mod_test_process_server('test_mod_process'))
