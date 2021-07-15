library(shiny)

ui <- fluidPage(
  uiOutput('tutu')
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    choice = 'test2'
  )

  output$tutu <- renderUI({
    selectInput('test', 'test',
                choices = c('None'='None', 'test1'='test1', 'test 2'='test2'),
                selected = rv$choice)
  })

}

shinyApp(ui, server)
