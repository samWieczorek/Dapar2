library(shiny)
library(R6)
library(Magellan)


source(file.path('.', 'Example_ProcessA.R'), local=TRUE)$value
source(file.path("../../R", "global.R"), local = TRUE)$value


actionBtnClass <- "btn-primary"

proc <- Example_ProcessA$new('processA', verbose = TRUE)

ui = fluidPage(
  tagList(
    proc$ui()
  )
  )

server = function(input, output){
  # Get a QFeatures dataset for example
  utils::data(Exp1_R25_prot, package='DAPARdata2')

  rv <- reactiveValues(
    dataIn = Exp1_R25_prot,
    dataOut = NULL
  )

  rv$dataOut <- proc$server(dataIn = reactive({Exp1_R25_prot}))
  }

shiny::shinyApp(ui, server)
