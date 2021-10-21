library(shiny)

setwd('~/GitHub/DaparToolshed/dev')

dirpath <- '../R'
for (l in list.files(path = dirpath, pattern = ".R"))
  source(file.path(dirpath, l), local=TRUE)$value


ui <- fluidPage(
  tagList(
    mod_popover_for_help_ui('popover'),
    selectInput("test", "", choices = seq_len(3))
  )
)



server <- function(input, output, session) {
  
  mod_popover_for_help_server("popover",
                              data = list(title = 'Test',
                                          content = "Help text.")
                              )
  }


shinyApp(ui=ui, server=server)
