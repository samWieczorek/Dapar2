#' @param id xxx
#' @export 
#' @importFrom shiny NS tagList 
#' @rdname descriptive-statistics
mod_ds_explorer_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("DS_sidebarPanel_tab")),
    uiOutput(ns("tabToShow"))
  )
}

#' @param id xxx
#' @param object An instance of the class `QFeatures`
#' @param which.assay A
#' @param digits xxx
#' 
#' 
#' @export
#' 
#' @return NA
#' @import DT
#' @importFrom tibble as_tibble
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment colData
#' 
#' @rdname descriptive-statistics
mod_ds_explorer_server <- function(id,
                                  object,
                                  which.assay,
                                  digits = reactive({3})
                                  ){ 
  
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    observe({
      req(object())
      stopifnot (inherits(object(), "QFeatures"))
    })
    
    
        output$DS_sidebarPanel_tab <- renderUI({
      
      .choices <- list( "Quantitative data" = "tabExprs",
                       "Proteins metadata" = "tabfData",
                       "Experimental design" = "tabpData")
      
      tagList(
        tags$div(
          tags$div( style="display:inline-block; vertical-align: middle; padding-right: 40px;",
                    radioButtons(ns("DS_TabsChoice"), 
                                 "Table to display",
                                 choices = .choices,
                                 inline = FALSE,
                                 selected = character(0))
          ),
          tags$div( style="display:inline-block; vertical-align: middle;",
                    uiOutput(ns("legend_ui"))
          )
          
        )
      )
    })
    
    
    
    
    mod_qMetadataLegend_server('legend')
    
    output$legend_ui <- renderUI({
      req(input$DS_TabsChoice == "tabExprs")
      mod_qMetadataLegend_ui(ns('legend'))
    })
    
    
    #----------------------------------------------
    output$tabToShow <- renderUI({
      req(input$DS_TabsChoice)
      switch(input$DS_TabsChoice,
             None = {return(NULL)},
             tabExprs = DT::DTOutput(ns("table")),
             tabfData = DT::DTOutput(ns("viewfData")),
             tabpData = DT::DTOutput(ns("viewDesign"))
             )
    })
    
    
    
    
    output$viewDesign <- DT::renderDT({
      req(object())
      
      data <- tibble::as_tibble(colData(object()))
      
      pal <- unique(RColorBrewer::brewer.pal(8, "Dark2"))
      
      dt <- DT::datatable(  data,
                            extensions = c('Scroller', 'Buttons'),
                            rownames=  FALSE,
                            options=list(initComplete = initComplete(),
                                         dom = 'Brtip',
                                         pageLength=10,
                                         orderClasses = TRUE,
                                         autoWidth=TRUE,
                                         deferRender = TRUE,
                                         bLengthChange = FALSE,
                                         scrollX = 200,
                                         scrollY = 500,
                                         scroller = TRUE,
                                         columnDefs = list(list(width='60px',targets= "_all"))
                            )) %>%
        DT::formatStyle(
          columns = colnames(data)[seq_len(2)],
          valueColumns = colnames(data)[2],
          backgroundColor = DT::styleEqual(unique(data$Condition), 
                                           pal[seq_len(length(unique(data$Condition)))])
        )
      
    })
    
    
    output$viewfData <- DT::renderDT({
      req(object())
      
      rdata <- rowData(object()[[which.assay()]])
      # Delete columns that are not one-dimensional
      rdata <- rdata[, - which(colnames(rdata) == 'adjacencyMatrix')]
      rdata <- rdata[, - which(colnames(rdata) == 'qMetadata')]
      
      dat <- DT::datatable(tibble::as_tibble(rdata),
                             rownames = TRUE,
                             extensions = c('Scroller', 'Buttons', 'FixedColumns'),
                             options=list(
                               initComplete = initComplete(),
                               dom = 'Bfrtip',
                               pageLength  = 10,
                               deferRender = TRUE,
                               bLengthChange = FALSE,
                               scrollX = 200,
                               scrollY = 600,
                               scroller = TRUE,
                               orderClasses = TRUE,
                               autoWidth = FALSE,
                               columns.searchable = FALSE,
                               fixedColumns = list(
                                 leftColumns = 1
                                 ),
                               columnDefs = list(
                                 list(
                                   columns.width = c("60px"),
                                   columnDefs.targets = c(list(0),list(1),list(2)))
                                 )
                               )
                           )

      if ('Significant' %in% colnames(rdata))
        dat <- dat %>%
        DT::formatStyle(columns = 'Significant',
                        target = 'row',
                        background = DT::styleEqual(1, 'lightblue'))
      
      return(dat)
    })
    
    
    output$table <- DT::renderDataTable(server=TRUE,{
      req(data())
      data <- object()[[which.assay()]]
      df <- cbind(keyId = rowData(data)[, idcol(data)],
                  round(assay(data), digits = digits()), 
                  qMetadata(data)
                  )
      mc <- qMetadata.def(typeDataset(data))
      colors <- as.list(setNames(mc$color, mc$node))
      
      DT::datatable( df,
                     extensions = c('Scroller'),
                     options = list(
                       initComplete = initComplete(),
                       displayLength = 20,
                       deferRender = TRUE,
                       bLengthChange = FALSE,
                       scrollX = 200,
                       scrollY = 600,
                       scroller = TRUE,
                       ordering = FALSE,
                       server = TRUE,
                       columnDefs = list(
                         list(
                           targets = c((( 2 + (ncol(df)-1)/2)):ncol(df)), 
                           visible = FALSE)
                         )
                       )
                     ) %>%
        formatStyle(
          colnames(df)[2:(1 + (ncol(df)-1)/2)],
          colnames(df)[((2 + (ncol(df)-1)/2)):ncol(df)],
          backgroundColor = styleEqual(names(colors), unlist(colors)),
          backgroundSize = '98% 48%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
      
    })
    

  
    #################
    # output$table <- DT::renderDT({
    #   req(object())
    #   df <- getDataForExprs()
    #   
    #   dt <- DT::datatable( df,
    #                        rownames=TRUE,
    #                        extensions = c('Scroller', 'Buttons', 'FixedColumns'),
    #                        options = list(
    #                          dom = 'Bfrtip',
    #                          initComplete = initComplete(),
    #                          displayLength = 20,
    #                          deferRender = TRUE,
    #                          bLengthChange = FALSE,
    #                          scrollX = 200,
    #                          scrollY = 600,
    #                          scroller = TRUE,
    #                          ordering=FALSE,
    #                          server = TRUE,
    #                          fixedColumns = list(leftColumns = 1),
    #                          columnDefs = list(list(targets = c(((ncol(df)/2)+1):ncol(df)), visible = FALSE)))) %>%
    #     DT::formatStyle(
    #       colnames(df)[seq_len(ncol(df)/2)],
    #       colnames(df)[seq((ncol(df)/2)+1, ncol(df))],
    #       #backgroundColor = DT::styleEqual(c("POV", "MEC"), c(rv.core$settings()$colorsTypeMV$POV, rv.core$settings()$colorsTypeMV$MEC)),
    #       backgroundColor = DT::styleEqual(c("POV", "MEC"), c("lightblue", "#E97D5E")), #orangeProstar)),
    #       backgroundSize = '98% 48%',
    #       backgroundRepeat = 'no-repeat',
    #       backgroundPosition = 'center'
    #     )
    #   
    #   dt
    # })
    
    
    
    getDataForExprs <- reactive({
      req(object())
      browser()
      test.table <- round(assay(object(), which.assay()), digits=10)
      test.table <- tibble::as_tibble(test.table)
      
      addon.table <- matrix(rep(NA,
                                ncol(test.table) * nrow(test.table)), 
                            nrow = nrow(test.table)
                            )
      addon.table <- tibble::as_tibble(addon.table)
      test.table <- cbind(test.table, addon.table)

      test.table
      
    })
    
    initComplete <- function(){
      
      return (DT::JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': 'darkgrey', 'color': 'black'});",
        "}"))
    }
    
    
  })
  
  
  
}
