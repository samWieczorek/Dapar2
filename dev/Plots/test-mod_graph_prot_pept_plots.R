library(visNetwork)
library(highcharter)


setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value


ui <- fluidPage(
  mod_graph_pept_prot_ui('plots_cc')
)


server <- function(input, output, session) {
  
  utils::data(Exp1_R25_pept, package='DAPARdata2')

  obj <- Exp1_R25_prot[[length(names(Exp1_R25_prot))]]
  
  mod_graph_pept_prot_server('plots_cc', 
                               obj = reactive({obj}),
                               cc = reactive({xxxx}),
                               matAdj = reactive({xxx})
  )
}


shinyApp(ui=ui, server=server)
