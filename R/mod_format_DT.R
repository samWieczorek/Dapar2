#' @title   mod_format_DT_ui and mod_format_DT_server
#' 
#' @description 
#' 
#' A shiny Module.
#'
#' @name mod_format_DT
#' 
#' 
#' @examples 
#' if(interactive()){
#'  library(QFeatures)
#'  library(shiny)
#'  library(DaparToolshed)
#'  data(ft)
#'  
#'   
#'  ui <- mod_format_DT_ui('dt')
#' 
#'  server <- function(input, output, session) {
#'   mod_format_DT_server(id = 'dt',
#'                        df = reactive({assay(ft,1)})
#'                        )
#'   }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
#' 
#' 
NULL


#' @param id shiny id
#' 
#' @importFrom shiny NS tagList 
#' @importFrom DT DTOutput
#' 
#' @export
#' @rdname mod_format_DT
#' 
mod_format_DT_ui <- function(id){
  ns <- NS(id)
  DT::DTOutput(ns("dt"))
}

#' @param id internal
#' @param df internal
#' @param rownames xxxx
#' @param dom xxx
#' @param style xxx
#' 
#' @export
#' 
#' @importFrom htmlwidgets JS    
#' @importFrom DT replaceData dataTableProxy renderDT datatable formatStyle styleEqual
#' @rdname mod_format_DT
mod_format_DT_server <- function(id,
                                 df,
                                 rownames = FALSE,
                                 dom = 'Bt',
                                 style = reactive({NULL})){
  
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    observe({
      req(df())
      DT::replaceData(proxy, df(), resetPaging = FALSE)  
    })
    
    proxy = DT::dataTableProxy(session$ns('dt'), session)
    
    output$dt <- DT::renderDT({
      req(df())
      isolate({
        if (is.null(style()) || length(style())==0){
          DT::datatable(df(), 
                        extensions = c('Scroller', 'Buttons'),
                        escape = FALSE,
                        rownames = rownames,
                        option = list(initComplete = initComplete(),
                                      dom = dom,
                                      server = FALSE,
                                      autoWidth = TRUE,
                                      columnDefs = list(list(width='150px', targets= "_all")),
                                      ordering = FALSE
                                      )
                        )
        } else {
          
          DT::datatable(df(), 
                        extensions = c('Scroller', 'Buttons'),
                        escape = FALSE,
                        rownames = rownames,
                        option = list(initComplete = initComplete(),
                                      dom = dom,
                                      server = FALSE,
                                      autoWidth = TRUE,
                                      columnDefs = list(list(width='150px',targets= "_all")),
                                      ordering = FALSE
                                      )
                        )  %>%
            DT::formatStyle(
              columns = style()$cols,
              valueColumns = style()$vals,
              backgroundColor = DT::styleEqual(style()$unique, 
                                               style()$pal)
            )
        }
      })
      
    })
    
    initComplete <- function(){
      
      return (htmlwidgets::JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': 'darkgrey', 'color': 'black'});",
        "}"))
    }
    
    
  })
  
}
