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
  
  
  bsCollapse(id = "collapseExample", 
             open = "",
             bsCollapsePanel(title = "Legend of colors",
                             uiOutput(ns('legend')),
                             style = ""
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
       
      output$legend <- renderUI({
        mc <- combine(qMetadata.def(level = 'peptide'),
                      qMetadata.def(level = 'protein')
        )
        
        tagList(
          lapply(1:nrow(mc), function(x){
            if (mc[x, 'color'] != 'white' || (mc[x, 'color'] == 'white' && !isTRUE(hide.white))) {
              tagList(
                tags$div(class="color-box",
                         style = paste0("display:inline-block; 
                                      vertical-align: middle;
                                      width:20px; height:20px;
                                      border:1px solid #000; 
                                      background-color: ", 
                                        mc[x, 'color'] , ";"),
                ),
                tags$p(style = paste0("display:inline-block; 
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
