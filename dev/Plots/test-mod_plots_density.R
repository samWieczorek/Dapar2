library(highcharter)
library(shiny)
library(SummarizedExperiment)


source(file.path("../../../R","mod_plots_density.R"), local=TRUE)$value
source(file.path("../../../R","mod_popover_for_help.R"), local=TRUE)$value

ui <- fluidPage(
  mod_plots_density_ui('plots_density')
)



server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- Exp1_R25_prot[[2]]
  conds <- colData(Exp1_R25_prot)[["Condition"]]
  legend <- colData(Exp1_R25_prot)[["Sample.name"]]

  mod_plots_density_server('plots_density',
                           obj = reactive({obj}),
                           conds = reactive({conds}),
                           legend = reactive({legend}),
                           base_palette = reactive({DaparToolshed::Example_Palette()})
  )
}


shinyApp(ui=ui, server=server)
