#' @param id xxx
#'
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom shiny NS tagList
#' 
#' @rdname q_metadata
#'
#' @export 
mod_qMetadataLegend_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    bsCollapse(id = "collapseExample", 
             open = "",
             bsCollapsePanel(title = "Legend of colors",
                             uiOutput(ns('legend')),
                             style = "info"
             )
  )
  )
  
}


#' @param id The 'id' of the shiny module.
#' @param hide.white A `logical(1)` that indicate if xxx
#' 
#' @export
#' 
#' @rdname q_metadata
mod_qMetadataLegend_server <- function(id, 
                                       hide.white = TRUE
                                       ){
  
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
       
      
      output$genericPlot <- renderPlot(plot(rnorm(100)))
      
      output$legend <- renderUI({
        mc <- combine(qMetadata.def(level = 'peptide'),
                      qMetadata.def(level = 'protein')
        )
        
        mc <- qMetadata.def(level = 'peptide')
        tagList(
          lapply(1:nrow(mc), function(x){
            cond <- mc[x, 'color'] == 'white' && hide.white
            cond <- cond || mc[x, 'color'] != 'white'
            if (cond) {
              tagList(
                tags$div(class="color-box",
                         style = paste0("display: inline-block; 
                                      vertical-align: middle;
                                      width: 20px; 
                                      height: 20px;
                                      border: 1px solid #000; 
                                      background-color: ", 
                                        mc[x, 'color'] , ";"),
                ),
                tags$p(style = paste0("display: inline-block; 
                                       vertical-align: middle;"),
                       mc[x, 'node']),
                br()
              )
            }
          })
        )
      })
    })
  
}
