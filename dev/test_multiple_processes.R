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

source(file.path('.', 'Example_ProcessA.R'), local=TRUE)$value
source(file.path('.', 'Example_ProcessB.R'), local=TRUE)$value
source(file.path('.', 'Example_Description.R'), local=TRUE)$value


#----------------------------------------------------------------------------





Pipeline <- R6Class(
  "Pipeline",
  public = list(
    id = NULL,
    ns = NULL,
    tmp.return = reactiveValues(),
    rv = reactiveValues(dataIn = NULL),
    child.process = list(
      #Example_Description = NULL,
      #Example_ProcessA = NULL
      #Filtering = NULL
      Normalization = NULL
      #Imputation = NULL
      #Example_ProcessB = NULL
    ),
    initialize = function(id){
      self$id <- id
      self$ns <- NS(id)
      self$child.process <- setNames(lapply(names(self$child.process),
                                            function(x){
                                              assign(x, base::get(paste0('Protein','_',x)))$new(self$ns(x))
                                            }),
                                     names(self$child.process)
      )
      # browser()
    },


    Launch_Servers = function(data){
      lapply(names(self$child.process), function(x){
        self$tmp.return[[x]] <- self$child.process[[x]]$server(dataIn = reactive({data()}))
      })
    },

    ui = function() {
      tagList(
        lapply(names(self$child.process), function(x){
          wellPanel(h3(x), self$child.process[[x]]$ui())
        })
      )
    },
    server = function(dataIn ) {

      self$Launch_Servers(data = reactive({dataIn()}))

      moduleServer(self$id, function(input, output, session) {

        output$show_ui <- renderUI({
          tagList(
            lapply(names(self$child.process), function(x){
              wellPanel(h3(x), self$child.process[[x]]$ui())
            })
          )
        })



      })
    }
  )
)

rv <- reactiveValues()
#Pipeline <- Pipeline$new('App')
process <- Protein_Normalization$new('Norm')
ui = fluidPage(
  tagList(
    actionButton('send', 'Send dataset'),
    actionButton('sendOldDataset', 'Send old dataset'),
    #Pipeline$ui()
    process$ui()
  )
)




server = function(session, input, output){
  # Get a QFeatures dataset for example

  rv <- reactiveValues(
    res = NULL
  )
  #rv$res <- Pipeline$server(dataIn = reactive({rv$dataIn}))
  rv$res <- process$server(dataIn = reactive({rv$dataIn}))
  utils::data(Exp1_R25_prot, package='DAPARdata2')


  # observe({
  #   req(rv$res())
  #  # print(metadata(rv$res()[[length(rv$res())]])$Params)
  # })

  observeEvent(input$send,{
    if (input$send%%2 != 0)
      rv$dataIn <- Exp1_R25_prot
    else
      rv$dataIn <- NULL
  })

  observeEvent(input$sendOldDataset,{
    if (input$sendOldDataset%%2 != 0){
      dataset <- addAssay(Exp1_R25_prot, Exp1_R25_prot[[length(Exp1_R25_prot)]], 'proteins_norm')
      metadata(dataset[[length(dataset)]])$Params <- list(
        method = 'GlobalQuantileAlignment',
        param2 = 'ABC',
        param3 = 36
      )
      rv$dataIn <- dataset
    }
    else
      rv$dataIn <- NULL
  })

}
shiny::shinyApp(ui, server)
