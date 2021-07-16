library(shiny)
library(highcharter)
library(SummarizedExperiment)


source(file.path("../../../R", "mod_plots_group_mv.R"), local=TRUE)$value
source(file.path('../../../R', 'mod_popover_for_help.R'), local=TRUE)$value
source(file.path('../../../R', 'global.R'), local=TRUE)$value


ui <- fluidPage(
  mod_plots_group_mv_ui('plots_group_mv')
)



server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- Exp1_R25_prot[[2]]

  conds <- SummarizedExperiment::colData(Exp1_R25_prot)

  mod_plots_group_mv_server('plots_group_mv',
                            obj = reactive({obj}),
                            conds = reactive({conds}),
                            base_palette = reactive({Dapar2::Example_Palette(conds, 
                                                                             Dapar2::Base_Palette(conditions = conds))})
  )

}


shinyApp(ui, server)
