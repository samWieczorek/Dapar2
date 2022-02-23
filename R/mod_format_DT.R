#' @title   mod_format_DT_ui and mod_format_DT_server
#' 
#' @description 
#' 
#' A shiny Module.
#'
#' @name mod_format_DT
#' 
NULL


#' @param id shiny id
#' 
#' @importFrom shiny NS tagList 
#' @importFrom DT renderDT DTOutput formatStyle %>% styleEqual
#' 
#' @export
#' @rdname mod_format_DT
#' 
mod_format_DT_ui <- function(id){
  ns <- NS(id)
  tagList(
    div(
      div( style="display:inline-block; vertical-align: middle; align: center;",
           DT::DTOutput(ns("dt"))
      )
    )
  )
}

#' @param id internal
#' @param table2show internal
#' @param rownames xxxx
#' @param dom xxx
#' @param style xxx
#' 
#' @export
#' 
#' @import DT
#' @importFrom htmlwidgets JS    
#' @importFrom DT replaceData dataTableProxy renderDT datatable formatStyle styleEqual
#' 
#' @return NA
#' @rdname mod_format_DT
mod_format_DT_server <- function(id,
                                 df,
                                 rownames = FALSE,
                                 dom = 'Bt',
                                 style = NULL){
  
  
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
                        rownames= rownames,
                        option=list(initComplete = initComplete(),
                                    dom = dom,
                                    server = FALSE,
                                    autoWidth=TRUE,
                                    columnDefs = list(list(width='150px',targets= "_all")),
                                    ordering = FALSE
                        )
          )
        } else {
          
          DT::datatable(df(), 
                        extensions = c('Scroller', 'Buttons'),
                        escape = FALSE,
                        rownames= rownames,
                        option=list(initComplete = initComplete(),
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
              backgroundColor = DT::styleEqual(style()$unique, style()$pal)
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
