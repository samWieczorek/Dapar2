library(shiny)

ui <- fluidPage(
  tagList(
    selectInput('select.1', 'Select 1', choice s= seq_len(4), width='100px'),
    selectInput('select.2', 'Select 2', choices = seq_len(4), width='100px'),
    selectInput('select.3', 'Select 3', choices = seq_len(4), width='100px'),
    selectInput('select.4', 'Select 4', choices = seq_len(4), width='100px'),
    uiOutput('show')

  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    widgets = list(
      select.1 = NULL,
      select.2 = NULL,
      select.3 = NULL,
      select.4 = NULL)
  )

observe({
  browser()
  lapply(seq_len(4), function(x){
    rv$widgets[[paste0('select.', x)]] <- input[[paste0('select.', x)]]
  })
})


  output$show <- renderUI({
    tagList(
      p(paste0(rv$widgets$select.1, ' - ', rv$widgets$select.2, ' - ', rv$widgets$select.3, ' - ', rv$widgets$select.4, ' - '))
    )
  })

}

shinyApp(ui, server)
