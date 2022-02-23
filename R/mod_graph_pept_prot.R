#' @title Connected components tools.
#'
#' @description 
#' 
#' A shiny Module.
#' 
#' @name connected_components
#' 
#' 
NULL

#' @param id Internal parameters for {shiny}.
#' 
#' @importFrom shiny NS tagList
#' @import visNetwork
#' @importFrom DT renderDT DTOutput formatStyle %>% styleEqual datatable
#' @importFrom shinyjs toggle hidden
#' 
#' @rdname connected_components
#' @export
mod_graph_pept_prot_ui <- function(id){
  ns <- NS(id)
  tagList(
    tabPanel("Graph",
             value = "graphTab",
             tabsetPanel(
               id = "graphsPanel",
               tabPanel("Settings",
                        uiOutput(ns("pepInfoUI"))
               ),
               tabPanel("CC one prot",
                        tagList(
                          shinyBS::bsCollapse(id = "collapseCCInfos",
                                              open = "",
                                              multiple = TRUE,
                                              shinyBS::bsCollapsePanel("One - One CC",
                                                                       fluidRow(
                                                                         column(width=4, DT::DTOutput(ns("OneOneDT"))),
                                                                         column(width=8, DT::DTOutput(ns("OneOneDTDetailed")))
                                                                       ),style = "info"),
                                              shinyBS::bsCollapsePanel("One - Multi CC",
                                                                       fluidRow(
                                                                         column(width=4, DT::DTOutput(ns("OneMultiDT"))),
                                                                         column(width=8, DT::DTOutput(ns("OneMultiDTDetailed")))
                                                                       ), style = "primary")
                          )
                        )
               ),
               tabPanel("CC multi prot",
                        tagList(
                          #uiOutput(ns("CCTooltip_UI")),
                          # highchartOutput(ns("jiji")))
                          selectInput(ns("searchCC"), 'Search for CC',
                                      choices = c('Tabular view' = 'tabular',
                                                  'Graphical view' = 'graphical'),
                                      width='150px'),
                          fluidRow(
                            column(width=6,tagList(
                              highchartOutput(ns("jiji")),
                              shinyjs::hidden( DTOutput(ns('CCMultiMulti')))
                            )),
                            column(width=6, tagList(
                              visNetworkOutput(ns("visNet_CC"), height='600px')))
                          ),
                          uiOutput(ns('CCDetailed'))
                        )
               )
             )
    )
  )
}


#' @param cc xxx
#' @param matAdj xxx
#' @param dataIn An object of class [QFeatures]
#' 
#' @rdname connected_components
#' @export
mod_graph_pept_prot_server <- function(id, cc, matAdj, dataIn){

  moduleServer(id, function(input, output, session){
    ns <- session$ns


    rv.cc <- reactiveValues(
      ## selected CC in global CC list (tab or plot)
      selectedCC = NULL,

      selectedNode = NULL,
      selectedNeighbors = NULL,
      selectedCCgraph = NULL,

      # when the user selects a node in the graph
      detailedselectedNode = list(
        sharedPepLabels = NULL,
        specPepLabels  = NULL,
        protLabels  = NULL)

    )



    output$pepInfoUI <- renderUI({
      dataIn()
      selectInput(ns('pepInfo'), "PepInfo", choices=colnames(Biobase::fData(dataIn())),
                  multiple=TRUE)
    })

    observeEvent(req(input$searchCC), {
      shinyjs::toggle('jiji', condition = input$searchCC=='graphical')
      shinyjs::toggle('CCMultiMulti', condition = input$searchCC=='tabular')
    })





    # select a point in the grpah
    observeEvent(input$click,{
      rv.cc$selectedNode <- input$click
    })


    # Get the id of selected neighbors in the graph
    observeEvent(input$visNet_CC_highlight_color_id,{
      rv.cc$selectedNeighbors <- input$visNet_CC_highlight_color_id
    })


    # select a CC in the summary table
    observeEvent(input$CCMultiMulti_rows_selected, {
      rv.cc$selectedCC <- input$CCMultiMulti_rows_selected
    })

    # select a CC in the jitter plot
    observeEvent(req(input$eventPointClicked), {
      this.index <- as.integer(strsplit(input$eventPointClicked, "_")[[1]][1])
      this.index+1
      rv.cc$selectedCC <- this.index+1
    })


    output$visNet_CC <- renderVisNetwork({
      req(rv.cc$selectedCC)
      local <- cc()[Get_CC_Multi2Any()]
      print(local[[rv.cc$selectedCC]])


      rv.cc$selectedCCgraph <- buildGraph(local[[rv.cc$selectedCC]], matAdj())

      display.CC.visNet(rv.cc$selectedCCgraph) %>%
        visEvents(click = paste0("function(nodes){
                Shiny.onInputChange('",ns("click"),"', nodes.nodes[0]);
                Shiny.onInputChange('",ns("node_selected"), "', nodes.nodes.length);
                ;}")
        ) %>%
        visOptions(highlightNearest = TRUE )
    })


    output$jiji <- renderHighchart({
      tooltip <- NULL

      isolate({
        local <-   cc()[Get_CC_Multi2Any()]
        n.prot <- unlist(lapply(local, function(x){length(x$proteins)}))
        n.pept <- unlist(lapply(local, function(x){length(x$peptides)}))
        df <- tibble(x=jitter(n.pept),
                     y = jitter(n.prot),
                     index = 1:length(local))

        if (!is.null( tooltip)){
          df <- cbind(df,Biobase::fData(dataIn())[ tooltip])
        }

        colnames(df) <- gsub(".", "_", colnames(df), fixed=TRUE)
        if (ncol(df) > 3){
          colnames(df)[4:ncol(df)] <-
            paste("tooltip_", colnames(df)[4:ncol(df)], sep="")
        }

        clickFun <-
          JS(paste0("function(event) {Shiny.onInputChange('",ns("eventPointClicked"),"', [this.index]+'_'+ [this.series.name]);}"))

        rv.core$tempplot$plotCC <-  DaparToolshed::plotJitter_rCharts(df, clickFunction=clickFun)

      })
      rv.core$tempplot$plotCC

    })



    output$CCMultiMulti <- DT::renderDT({
      Get_CC_Multi2Any()
      df <- do.call(rbind,lapply(cc()[Get_CC_Multi2Any()],
                                 function(x){
                                   data.frame(rbind(x),
                                              nPep = length(x$peptides),
                                              nProt = length(x$proteins))}))
      df <- cbind(df,id = 1:nrow(df))
      df <- df[c('id', 'nProt', 'nPep', 'proteins', 'peptides')]

      dat <- DT::datatable(df,
                           selection = 'single',
                           rownames=FALSE,
                           extensions = c('Scroller', 'Buttons'),
                           options=list(initComplete = initComplete(),
                                        dom='Bfrtip',
                                        deferRender = TRUE,
                                        bLengthChange = FALSE,
                                        scrollX = 400,
                                        scrollY = 400,
                                        displayLength = 10,
                                        scroller = TRUE,
                                        orderClasses = TRUE,
                                        autoWidth=TRUE,
                                        columns.searchable = FALSE,
                                        columnDefs = list(list(columns.width=c("60px"),
                                                               columnDefs.targets=c(list(0),list(1),list(2))))))

      return(dat)
    })



    observeEvent(c(rv.cc$selectedNeighbors,input$node_selected,rv.cc$selectedCCgraph), {

      local <-   cc()[Get_CC_Multi2Any()]
      rv.cc$selectedNeighbors

      nodes <- rv.cc$selectedCCgraph$nodes

      if(!is.null(input$node_selected) && input$node_selected == 1){
        sharedPepIndices <- intersect(rv.cc$selectedNeighbors, which(nodes[,'group'] == "shared.peptide"))
        specPepIndices <- intersect(rv.cc$selectedNeighbors, which(nodes[,'group'] == "spec.peptide"))
        protIndices <- intersect(rv.cc$selectedNeighbors,which(nodes[,'group'] == "protein"))

      } else {
        sharedPepIndices <- which(nodes[,'group'] == "shared.peptide")
        specPepIndices <- which(nodes[,'group'] == "spec.peptide")
        protIndices <- which(nodes[,'group'] == "protein")
      }
      rv.cc$detailedselectedNode$sharedPepLabels <- as.numeric(nodes[sharedPepIndices, 'label'])
      rv.cc$detailedselectedNode$specPepLabels <-  as.numeric(nodes[specPepIndices, 'label'])
      rv.cc$detailedselectedNode$protLabels <-  as.numeric(nodes[protIndices, 'label'])

    })


    output$CCDetailed <- renderUI({
      req(rv.cc$detailedselectedNode)
      req(rv.cc$selectedCC)

      tagList(
        h4("Proteins"),
        DT::DTOutput(ns('CCDetailedProt')),
        h4("Specific peptides"),
        DT::DTOutput(ns('CCDetailedSpecPep')),
        h4("Shared peptides"),
        DT::DTOutput(ns('CCDetailedSharedPep'))
      )
    })

    output$CCDetailedProt<- DT::renderDT({
      req(rv.cc$selectedCC)
      rv.cc$detailedselectedNode
      if(is.null(rv.cc$detailedselectedNode$protLabels)){return(NULL)}

      df <- data.frame(proteinId = unlist(rv.cc$detailedselectedNode$protLabels)
                       #other = rep(NA,length(rv.cc$detailedselectedNode$protLabels))
      )
      dt <- datatable( df,
                       extensions = c('Scroller'),
                       options = list(initComplete = initComplete(),
                                      dom='rt',
                                      blengthChange = FALSE,
                                      ordering=FALSE,
                                      scrollX = 400,
                                      scrollY = 100,
                                      displayLength = 10,
                                      scroller = TRUE,
                                      header=FALSE,
                                      server = FALSE))
      dt
    })



    #######

    output$CCDetailedSharedPep <- DT::renderDT({
      rv.cc$detailedselectedNode
      input$pepInfo

      if(is.null((rv.cc$detailedselectedNode$sharedPepLabels))){return(NULL)}


      ind <- 1:ncol(dataIn())
      data <- getDataForExprs(dataIn())
      pepLine <- rv.cc$detailedselectedNode$sharedPepLabels
      indices <- unlist(lapply(pepLine, function(x){which(rownames(data)==x)}))
      data <- data[indices,c(ind, (ind + ncol(data)/2))]

      if(!is.null(input$pepInfo))
      {
        data <- cbind(data, Biobase::fData(dataIn())[pepLine,input$pepInfo])
        colnames(data)[(1+ncol(data)-length(input$pepInfo)):ncol(data)] <- input$pepInfo
      }

      offset <- length(input$pepInfo)
      dt <- datatable( data,
                       extensions = c('Scroller'),
                       options = list(initComplete = initComplete(),
                                      dom='rt',
                                      blengthChange = FALSE,
                                      ordering=FALSE,
                                      scrollX = 400,
                                      scrollY = 150,
                                      displayLength = 10,
                                      scroller = TRUE,
                                      header=FALSE,
                                      server = FALSE,
                                      columnDefs = list(list(targets = c((((ncol(data)-offset)/2)+1):(ncol(data)-offset)), visible = FALSE))
                       )) %>%
        DT::formatStyle(
          colnames(data)[1:((ncol(data)-offset)/2)],
          colnames(data)[(((ncol(data)-offset)/2)+1):(ncol(data)-offset)],
          backgroundColor = DT::styleEqual(c("POV", "MEC"), c(rv.prostar$settings()$colorsTypeMV$POV, rv.prostar$settings()$colorsTypeMV$MEC)))

      dt
    })





    #####-----------
    output$CCDetailedSpecPep <- DT::renderDT({
      rv.cc$detailedselectedNode
      input$pepInfo
      if(is.null((rv.cc$detailedselectedNode$specPepLabels))){return(NULL)}

      ind <- 1:ncol(dataIn())
      data <- getDataForExprs(dataIn())
      pepLine <-  rv.cc$detailedselectedNode$specPepLabels
      indices <- unlist(lapply(pepLine, function(x){which(rownames(data)==x)}))
      data <- data[indices,c(ind, (ind + ncol(data)/2))]

      if(!is.null(input$pepInfo))
      {
        data <- cbind(data, Biobase::fData(dataIn())[pepLine,input$pepInfo])
        colnames(data)[(1+ncol(data)-length(input$pepInfo)):ncol(data)] <- input$pepInfo
      }

      offset <- length(input$pepInfo)


      dt <- datatable( data,
                       extensions = c('Scroller'),
                       options = list(initComplete = initComplete(),
                                      dom='rt',
                                      blengthChange = FALSE,
                                      ordering=FALSE,
                                      scrollX = 400,
                                      scrollY = 100,
                                      displayLength = 10,
                                      scroller = TRUE,
                                      header=FALSE,
                                      server = FALSE,
                                      columnDefs = list(list(targets = c((((ncol(data)-offset)/2)+1):(ncol(data)-offset)), visible = FALSE))
                       )) %>%
        DT::formatStyle(
          colnames(data)[1:((ncol(data)-offset)/2)],
          colnames(data)[(((ncol(data)-offset)/2)+1):(ncol(data)-offset)],
          backgroundColor = DT::styleEqual(c("POV", "MEC"), c(rv.prostar$settings()$colorsTypeMV$POV, rv.prostar$settings()$colorsTypeMV$MEC)))

      dt
    })






    Get_CC_One2One <- reactive({
      cc()
      ll.prot <- lapply(cc(), function(x){length(x$proteins)})
      ll.pept <- lapply(cc(), function(x){length(x$peptides)})
      ll.prot.one2one <- intersect(which(ll.prot == 1),which(ll.pept == 1))
      ll.prot.one2one
    })

    Get_CC_One2multi <- reactive({
      cc()
      ll.prot <- lapply(cc(), function(x){length(x$proteins)})
      ll.pept <- lapply(cc(), function(x){length(x$peptides)})
      ll.prot.one2multi <- intersect(which(ll.prot == 1),which(ll.pept > 1))
      ll.prot.one2multi
    })

    Get_CC_Multi2Any <- reactive({
      cc()
      ll.prot <- lapply(cc(), function(x){length(x$proteins)})
      ll.pept <- lapply(cc(), function(x){length(x$peptides)})
      ll.prot.multi2any <- which(ll.prot > 1)
      ll.prot.multi2any
    })



    BuildOne2OneTab <- reactive({
      cc()
      table <- do.call(rbind,lapply(cc()[Get_CC_One2One()],function(x){data.frame(rbind(x))}))
      table
    })

    BuildOne2MultiTab <- reactive({
      cc()
      table <- do.call(rbind,lapply(cc()[Get_CC_One2multi()],function(x){data.frame(rbind(x), nPep = length(x$peptides))}))
      table <- table[c('proteins', 'nPep', 'peptides')]
      table
    })


    BuildMulti2AnyTab <- reactive({
      cc()
      table <- do.call(rbind,lapply(cc()[Get_CC_Multi2Any()],function(x){data.frame(rbind(x), nPep = length(x$peptides))}))
      table <- table[c('proteins', 'nPep', 'peptides')]

      table
    })





    output$OneMultiDT <- DT::renderDT({
      req(cc())

      dat <- DT::datatable(BuildOne2MultiTab(),
                           selection = 'single',
                           rownames=FALSE,
                           extensions = c('Scroller', 'Buttons'),
                           options=list(initComplete = initComplete(),
                                        dom='Bfrtip',
                                        deferRender = TRUE,
                                        bLengthChange = TRUE,
                                        displayLength = 10,
                                        scrollX = 400,
                                        scrollY = 400,
                                        scroller = TRUE,
                                        orderClasses = TRUE,
                                        autoWidth = FALSE,
                                        columns.searchable = FALSE,
                                        columnDefs = list(list(columns.width=c("60px"),
                                                               columnDefs.targets=c(list(0),list(1),list(2))))))

      return(dat)
    })



    output$OneMultiDTDetailed <- DT::renderDT({
      input$pepInfo
      req(input$OneMultiDT_rows_selected)

      line <- input$OneMultiDT_rows_selected

      ind <- 1:ncol(dataIn())
      data <- getDataForExprs(dataIn())
      pepLine <- as.numeric(unlist(BuildOne2MultiTab()[line,"peptides"]))

      indices <- unlist(lapply(pepLine, function(x){which(rownames(data)==x)}))

      data <- data[indices,c(ind, (ind + ncol(data)/2))]

      if(!is.null(input$pepInfo))
      {
        data <- cbind(data, Biobase::fData(dataIn())[pepLine,input$pepInfo])
        colnames(data)[(1+ncol(data)-length(input$pepInfo)):ncol(data)] <- input$pepInfo
      }

      offset <- length(input$pepInfo)


      dt <- datatable( data,
                       extensions = c('Scroller', 'Buttons'),
                       options = list(initComplete = initComplete(),
                                      dom='Bfrtip',
                                      pageLength = 10,
                                      blengthChange = FALSE,
                                      displayLength = 10,
                                      ordering=FALSE,
                                      header=FALSE,
                                      server = FALSE,
                                      columnDefs = list(list(targets = c((((ncol(data)-offset)/2)+1):(ncol(data)-offset)), visible = FALSE))
                       )) %>%
        DT::formatStyle(
          colnames(data)[1:((ncol(data)-offset)/2)],
          colnames(data)[(((ncol(data)-offset)/2)+1):(ncol(data)-offset)],
          backgroundColor = DT::styleEqual(c("POV", "MEC"), c(rv.prostar$settings()$colorsTypeMV$POV, rv.prostar$settings()$colorsTypeMV$MEC)))

      dt
    })




    output$OneOneDT <- DT::renderDT({
      req(cc())

      dat <- DT::datatable(BuildOne2OneTab(),
                           selection = 'single',
                           rownames=FALSE,
                           extensions = c('Scroller', 'Buttons'),
                           options=list(initComplete = initComplete(),
                                        dom='Bfrtip',
                                        deferRender = TRUE,
                                        bLengthChange = FALSE,
                                        scrollX = 400,
                                        scrollY = 200,
                                        scroller = TRUE,
                                        orderClasses = TRUE,
                                        autoWidth = FALSE,
                                        columns.searchable = FALSE,
                                        columnDefs = list(list(columns.width=c("60px"),
                                                               columnDefs.targets=c(list(0),list(1),list(2))))))

      return(dat)
    })



    output$OneOneDTDetailed <- DT::renderDT({
      req(cc())
      req(input$OneOneDT_rows_selected)
      input$pepInfo

      line <- input$OneOneDT_rows_selected

      ind <- 1:ncol(dataIn())
      data <- getDataForExprs(dataIn())
      pepLine <- as.numeric(BuildOne2OneTab()[line,2])
      indices <- unlist(lapply(pepLine, function(x){which(rownames(data)==x)}))
      data <- data[indices,c(ind, (ind + ncol(data)/2))]

      if(!is.null(input$pepInfo))
      {
        data <- cbind(data, Biobase::fData(dataIn())[pepLine,input$pepInfo])
        colnames(data)[(1+ncol(data)-length(input$pepInfo)):ncol(data)] <- input$pepInfo
      }

      offset <- length(input$pepInfo)

      dt <- datatable( data,
                       extensions = c('Scroller', 'Buttons'),
                       options = list(initComplete = initComplete(),
                                      dom='Bfrtip',
                                      blengthChange = FALSE,
                                      pageLength = 10,
                                      displayLength = 10,
                                      ordering=FALSE,
                                      header=FALSE,
                                      server = FALSE,
                                      columnDefs = list(list(targets = c((((ncol(data)-offset)/2)+1):(ncol(data)-offset)), visible = FALSE))
                       )) %>%
        DT::formatStyle(
          colnames(data)[1:((ncol(data)-offset)/2)],
          colnames(data)[(((ncol(data)-offset)/2)+1):(ncol(data)-offset)],
          backgroundColor = DT::styleEqual(c("POV", "MEC"), c(rv.prostar$settings()$colorsTypeMV$POV, rv.prostar$settings()$colorsTypeMV$MEC)))

      dt
    })

  })


}

## To be copied in the UI
# mod_graph_pept_prot_ui("graph_pept_prot_ui_1")

## To be copied in the server
# callModule(mod_graph_pept_prot_server, "graph_pept_prot_ui_1")
