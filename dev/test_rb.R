library(shiny)
library(base64enc)

ui <- fluidPage(
  tagList(
    tags$style(HTML("
            .shiny-input-radiogroup label {
                display: inline-block;
                text-align: center;
                margin-bottom: 20px;
            }
            .shiny-input-radiogroup label input[type='radio'] {
                display: inline-block;
                margin: 3em auto;
                margin-left: 10px;
            }
        ")),
  uiOutput('show_rb'),
  verbatimTextOutput(outputId = "res1"),
)
)

server <- function(input, output, session) {
  output$res1 <- renderPrint(input$rb)

  output$show_rb <- renderUI({

    ll <- c('mod_plots_var_dist', 'mod_plots_var_dist', 'mod_plots_var_dist')
    radioButtons("rb", "",
                 inline = T,
                 choiceNames = lapply(ll, function(x){
                   img(src = dataURI(file=system.file('images', paste0(x, '.png'), package="MSPipelines"), mime="image/png"),
                       width='30px')
                 }),
                 choiceValues = list(
                   "icon", "html", "text"
                 ))
  })
}

# Run the application
shinyApp(ui = ui, server = server)






#
#
#
#
# ui <- fluidPage(
#    radioButtons("rb", "",
#                 inline = T,
#                choiceNames = list(
#                  img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                      width='30px'),
#                  img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                      width='30px'),
#                  img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                      width='30px')
#                 ),
#                choiceValues = list(
#                  "icon", "html", "text"
#                 )),
#    textOutput("txt")
#  )
#
#  server <- function(input, output) {
#   output$txt <- renderText({
#      paste("You chose", input$rb)
#    })
#  }
#
#  shinyApp(ui, server)
#
#
#
#
#
#
#
#
#  library(shiny)
#  library(shinyWidgets)
#
#  ui <- fluidPage(
#    tags$h1("Pretty radio buttons"),
#    br(),
#
#    fluidRow(
#      tags$style(HTML("
#             .shiny-input-radiogroup label {
#                 display: inline-block;
#                 text-align: center;
#             }
#             .shiny-input-radiogroup label input[type='radio'] {
#                 display: block;
#                 margin: 2em auto;
#             }
#
#         ")),
#      column(
#        width = 4,
#        prettyRadioButtons(
#          inputId = "radio1",
#          label = "Click me!",
#          inline = T,
#          choiceNames = list(
#            img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                width='30px'),
#            img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                width='30px'),
#            img(src = base64enc::dataURI(file=system.file('images', 'mod_plots_var_dist.png', package="MSPipelines"), mime="image/png"),
#                width='30px')
#          ),
#          choiceValues = list(
#            "icon", "html", "text"
#          )
#        ),
#        verbatimTextOutput(outputId = "res1"),
#
#      )
#
#    )
#
#  )
#
#  server <- function(input, output, session) {
#
#    output$res1 <- renderPrint(input$radio1)
#
#
#  }
#
#  if (interactive())
#    shinyApp(ui, server)
#
