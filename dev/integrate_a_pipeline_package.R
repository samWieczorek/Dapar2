library(MSPipelines)
library(shiny)


## Example of Protein pipeline from MSPipelines package
pipeline <- Protein$new('protPipeline')

ui = fluidPage(
  tagList(
    actionButton('send', 'Send dataset'),
    pipeline$ui()
  )
)

server = function(input, output, session){

  # Get a QFeatures dataset for example
  basename(f <- msdata::quant(pattern = "cptac", full.names = TRUE))
  i <- grep("Intensity\\.", names(read.delim(f)))
  cptac <- QFeatures::readQFeatures(f, ecol = i, sep = "\t", name = "peptides", fnames = "Sequence")

  rv <- reactiveValues(
    res = NULL,
    dataIn = NULL
  )

  rv$res <- pipeline$server(dataIn = reactive({rv$dataIn}))

  GetResult <- reactive({
    triggerValues <- unlist(rv$res()$trigger)
    if (sum(triggerValues)>0){ # Init of core engine
      processHasChanged <- pipeline$config$steps[which(max(triggerValues)==triggerValues)]
      ind.processHasChanged <- which(pipeline$config$steps==processHasChanged)
      newValue <- rv$res()$value
      print(names(newValue))
    }
  })


  observe({ GetResult() })


  observeEvent(input$send,{
    if (input$send%%2 != 0)
      rv$dataIn <- cptac
    else
      rv$dataIn <- NULL
  })
}

shiny::shinyApp(ui, server)
