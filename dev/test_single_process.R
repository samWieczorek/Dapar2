library(shiny)
library(R6)
library(tibble)
library(Magellan)
library(MSPipelines)
library(highcharter)

options(shiny.fullstacktrace = T)

#------------------------ Class TimelineDraw --------------------------------------
source(file.path('../R', 'global.R'), local=TRUE)$value
source(file.path("../R", "mod_popover_for_help.R"), local = TRUE)$value
source(file.path("../R", "mod_format_DT.R"), local = TRUE)$value
source(file.path("../R", "Protein_Normalization.R"), local = TRUE)$value

# source(file.path('.', 'Example_ProcessA.R'), local=TRUE)$value
# source(file.path('.', 'Example_ProcessB.R'), local=TRUE)$value
# source(file.path('.', 'Example_Description.R'), local=TRUE)$value


rv <- reactiveValues()

process <- Protein_Normalization$new('Norm')

ui = fluidPage(
  tagList(
    actionButton('send', 'Send dataset'),
    actionButton('sendOldDataset', 'Send old dataset'),
    actionButton('updateStatus', 'Update status'),
    process$ui()
  )
)




server = function(session, input, output){
  # Get a QFeatures dataset for example

  rv <- reactiveValues(
    res = NULL
  )
  utils::data(Exp1_R25_prot, package='DAPARdata2')

  #rv$res <- process$server(dataIn = reactive({rv$dataIn}))
  rv$res <- process$server(dataIn = reactive({Exp1_R25_prot}))

  observe({
    req(rv$res()$value)
    print(metadata((rv$res()$value)[[length(rv$res())]])$Params)
  })

  observeEvent(input$updateStatus, {
    process$rv$status <- c(1, 1, 1)
  })

  observeEvent(input$send,{
    if (input$send%%2 != 0)
      rv$dataIn <- Exp1_R25_prot
    else
      rv$dataIn <- NULL
  })

  observeEvent(input$sendOldDataset,{
    if (input$sendOldDataset%%2 != 0){
      dataset <- addAssay(Exp1_R25_prot, Exp1_R25_prot[[length(Exp1_R25_prot)]], 'proteins_norm')
      metadata(dataset[[length(dataset)]])$status <- c(1, 1, 1)
      metadata(dataset[[length(dataset)]])$Params$Normalize <- list(
        method = "QuantileCentering",
        type = "within conditions",
        quantile = 0.15,
        typeSelect = "ProteinList",
        randSelect = "",
        colSelect = NULL,
        rand.indices = "",
        col.indices = NULL,
        listSelect = c("CON__P04264", "CON__P07477", "CON__P13645"),
        list.indices = c(4, 5, 6)
        )
      rv$dataIn <- dataset
    }
    else
      rv$dataIn <- NULL
  })

  }
shiny::shinyApp(ui, server)
