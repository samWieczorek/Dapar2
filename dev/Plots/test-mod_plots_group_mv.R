library(shiny)
library(highcharter)
library(SummarizedExperiment)


setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value


ui <- fluidPage(
  mod_plots_group_mv_ui('plots_group_mv')
)



server <- function(input, output, session) {

  utils::data(Exp1_R25_prot, package='DAPARdata2')

  obj <- Exp1_R25_prot[[2]]

  samplesData <- SummarizedExperiment::colData(Exp1_R25_prot)
  conds <- samplesData$Condition

  mod_plots_group_mv_server('plots_group_mv',
                            obj = reactive({obj}),
                            conds = reactive({conds}),
                            base_palette = reactive({DaparToolshed::Example_Palette(conds, 
                                                                                    DaparToolshed::Base_Palette(conditions = conds)
                                                                                    )
                              })
                            )

}


shinyApp(ui, server)
