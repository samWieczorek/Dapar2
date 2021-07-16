library(shiny)
library(highcharter)
library(SummarizedExperiment)


source(file.path('../../R', 'mod_plots_tracking.R'), local=TRUE)$value
source(file.path("../../R","mod_plots_intensity.R"), local=TRUE)$value
source(file.path("../../R","mod_popover_for_help.R"), local=TRUE)$value


ui <- fluidPage(
  uiOutput('check'),
  mod_plots_tracking_ui('master_tracking'),
  mod_plots_intensity_ui('plots_boxplots')
)



server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')
  keyId <- metadata(Exp1_R25_prot)[['keyId']]

  r <- reactiveValues(
    master = NULL,
    widgets = list(
      sync = FALSE
    )
  )

  metadata <- metadata(Exp1_R25_prot)
  conds <- colData(Exp1_R25_prot)[['Condition']]
  obj <- Exp1_R25_prot[[2]]
  SummarizedExperiment::rowData(obj) <- cbind(SummarizedExperiment::rowData(obj),
                                              ProtOfInterest=rep(0,nrow(obj)))
  SummarizedExperiment::rowData(obj)$ProtOfInterest[10:20] <- 1


  r$master <- mod_plots_tracking_server('master_tracking',
                                        obj = reactive({obj}),
                                        keyId=reactive({keyId}),
                                        params=reactive({NULL}),
                                        reset=reactive({FALSE}),
                                        slave = reactive({FALSE})
  )

  mod_plots_intensity_server('plots_boxplots',
                             dataIn = reactive({obj}),
                             meta = reactive({metadata}),
                             conds = reactive({conds}),
                             base_palette = reactive({Dapar2::Example_Palette(conds, Dapar2::Base_Palette(conditions = conds))}),
                             params = reactive({if(r$widgets$sync) r$master() else NULL }),
                             reset = reactive({FALSE}),
                             slave = reactive({r$widgets$sync})
  )


  observeEvent(input$sync, { r$widgets$sync <- input$sync})

  output$check <- renderUI({

    checkboxInput('sync', 'sync Slave with Master', value=r$widgets$sync)
  })


  }


shinyApp(ui, server)
