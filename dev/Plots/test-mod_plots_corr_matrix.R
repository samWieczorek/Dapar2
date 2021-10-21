library(shiny)
library(highcharter)
library(SummarizedExperiment)

setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value


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
