library(shiny)
library(highcharter)
library(SummarizedExperiment)

source(file.path("../../../R","mod_plots_corr_matrix.R"), local=TRUE)$value
source(file.path("../../../R","mod_popover_for_help.R"), local=TRUE)$value
source(file.path("../../../R","global.R"), local=TRUE)$value



ui <- fluidPage(
  mod_plots_corr_matrix_ui('plots_corr_matrix')
)


server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- Exp1_R25_prot[[length(names(Exp1_R25_prot))]]
  names <- gsub('Intensity_','',colnames(assay(Exp1_R25_prot)))

  mod_plots_corr_matrix_server('plots_corr_matrix',
             obj = reactive({obj}),
             names = reactive({NULL}),
             gradientRate = reactive({NULL})
             )
}


shinyApp(ui=ui, server=server)
