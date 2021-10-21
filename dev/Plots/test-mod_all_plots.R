library(shiny)
library(highcharter)

setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value



ui <- fluidPage(
  mod_all_plots_ui('plots')
)


server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- QFeatures::addAssay(Exp1_R25_prot, (QFeatures::filterNA(Exp1_R25_prot,i=2))[[2]], "original_log_NAfiltered")

  mod_all_plots_server('plots',
                       dataIn = reactive({obj})
                       )
}


shinyApp(ui, server)
