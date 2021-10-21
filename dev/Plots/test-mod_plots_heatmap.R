library(shiny)
library(SummarizedExperiment)

setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value


ui <- fluidPage(
  mod_plots_heatmap_ui('plots_heatmap')
)


server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- Exp1_R25_prot[[2]]
  conds <- colData(Exp1_R25_prot)[['Condition']]

  mod_plots_heatmap_server('plots_heatmap',
                           obj = reactive({obj}),
                           conds = reactive({conds})
  )

}


shinyApp(ui, server)
