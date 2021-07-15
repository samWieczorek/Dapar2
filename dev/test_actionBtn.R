library(shiny)

ui <- fluidPage(
  actionButton(inputId = "A.button",
               label = img(src=system.file("images", "desc_pca.png", package="MSPipelines"), width="30", height="30"),
               style = "width: 50px; height: 50px;
background: 'desc_pca.png';  background-size: cover; background-position: center;"),

     div( style="display:inline-block; vertical-align: middle; padding: 7px",
         tags$button(
           id = "btn_3",
           class = "btn action-button",
           img(src = system.file("extdata/images", "desc_pca.png", package="MSPipelines"),
               height = "30px")
         )
    )
)

server <- function(input, output, session) {


}

shinyApp(ui, server)
