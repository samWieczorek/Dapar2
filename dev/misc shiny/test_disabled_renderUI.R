library(shiny)
library(shinyjs)
ui <- fluidPage(
  useShinyjs(),
  uiOutput("id1"),
  actionButton("id2","No render")
)

server <- function(input, output, session) {

  output$id1=renderUI({
    actionButton("id0","renderUI")
  })

delay(1, shinyjs::disable("id0"))
  shinyjs::disable("id2")

}

shinyApp(ui, server)
