# Module UI
  
#' @title   mod_legend_colored_exprs_ui and mod_legend_colored_exprs_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname descriptive_statistics_plots
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
#' @return NA
#' 
mod_plots_legend_colored_exprs_ui <- function(id){
  ns <- NS(id)
  tagList(
    tags$p(tags$b("Legend of colors")),
    
    fluidRow(
      column(width=2, HTML(paste0("<div style=\"width:50px;height:20px;border:0px solid #000; background-color: ",
                                  'lightblue',";\"></div>")) ),
      column(width=10, tags$p("Partially Observed Value"))
    ),
    
    fluidRow(
      column(width=2,HTML(paste0("<div style=\"width:50px;height:20px;border:0px solid #000; background-color: ",
                                 #orangeProstar,";\"></div>"))),
                                 "#E97D5E",";\"></div>"))),
      column(width=10, tags$p("Missing in Entire Condition"))
    )
  )
}
    
# Module Server
    
#' @rdname descriptive_statistics_plots
#' @export
#' @keywords internal
#' @return NA
    
mod_plots_legend_colored_exprs_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
  })
  
}
    
## To be copied in the UI
# mod_plots_legend_colored_exprs_ui("legend_colored_exprs_ui_1")
    
## To be copied in the server
# callModule(mod_plots_legend_colored_exprs_server, "legend_colored_exprs_ui_1")
 
