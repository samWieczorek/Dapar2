library(shiny)
library(R6)
library(tibble)
library(Magellan)
library(MSPipelines)

verbose <- T

options(shiny.fullstacktrace = T)
options(shiny.reactlog=TRUE)

#------------------------ Class TimelineDraw -----------------------------------
# source(file.path('../R', 'global.R'), local=TRUE)$value
# source(file.path("../R", "mod_popover_for_help.R"), local = TRUE)$value
# source(file.path("../R", "mod_format_DT.R"), local = TRUE)$value
# source(file.path("../R", "Protein_Filtering.R"), local = TRUE)$value
# source(file.path('../R', 'Protein_Description.R'), local=TRUE)$value
# source(file.path('../R', 'Protein_Normalization.R'), local=TRUE)$value
# source(file.path('../R', 'Protein_Imputation.R'), local=TRUE)$value
# source(file.path('../R', 'Protein_HypothesisTest.R'), local=TRUE)$value
# source(file.path('../R', 'Protein.R'), local=TRUE)$value


## Main app
pipeline <- Protein$new('App')

ui = fluidPage(
  tagList(
    div( style="display:inline-block; vertical-align: middle; padding: 7px",
         actionButton('send', 'Send dataset')
    ),
    div( style="display:inline-block; vertical-align: middle; padding: 7px",
         actionButton('showPlotModal', 'Show plots')
    ),
    div( style="display:inline-block; vertical-align: middle; padding: 7px",
         mod_bsmodal_ui('plotModal')
    ),

  pipeline$ui()
  )
)

server = function(input, output, session){

  # Get a QFeatures dataset for example
  utils::data(Exp1_R25_prot, package='DAPARdata2')

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
      rv$dataIn <- Exp1_R25_prot
    else
      rv$dataIn <- NULL
  })

  MSPipelines::mod_all_plots_server('plots',
                                    dataIn = reactive({rv$dataIn})
  )
  mod_bsmodal_server('plotModal',
                     title = 'Plots',
                     width="75%", # en px ou % de largeur
                     uiContent = MSPipelines::mod_all_plots_ui('plots')
  )
}

shiny::shinyApp(ui, server)
