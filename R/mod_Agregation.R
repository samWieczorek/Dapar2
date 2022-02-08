
#' @title Module agregation
#' 
#' @param id
#' 
#' @rdname mod_Agregation
#' 
#' @export
#' 
mod_Agregation_ui <- function(id){
  ns <- NS(id)
}


#' @param id xxx
#' @param dataIn xxx
#' @param steps.enabled xxx
#' @param remoteReset xxx
#' @param steps.status xxx
#' @param current.pos xxx
#' @param verbose xxx
#' 
#' 
#' 
#' 
#' @rdname mod_Agregation
#' @importFrom shinyjs toggle hidden
#' 
#' @export
#' 
mod_Agregation_server <- function(id,
                                  dataIn = reactive({NULL}),
                                  steps.enabled = reactive({NULL}),
                                  remoteReset = reactive({FALSE}),
                                  steps.status = reactive({NULL}),
                                  current.pos = reactive({1}),
                                  verbose = FALSE
){
  
  
  # This list contains the basic configuration of the process
  config <- list(
    # Define the type of module
    mode = 'process',
    
    # List of all steps of the process
    steps = c('Description', 'Agregation', 'Save'),
    
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, TRUE, TRUE),
    
    path_to_md_dir =  system.file('md/', package='DaparToolshed')
  )
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    Agregation_ProteinId = NULL,
    Agregation_useOfShared = NULL,
    Agregation_consider = NULL,
    Agregation_nTopn = NULL,
    Agregation_operator = NULL,
    nbPeptides = 0,
    Agregation_filterProtAfterAgregation = FALSE,
    Agregation_barplotType = 'all',
    AddMetadata_columnsForProteinDatasetBox = NULL
  )
  
  
  
  
  
  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(str2expression(Get_Worflow_Core_Code(
      w.names = names(widgets.default.values)
    )))
    
    rv.custom <- reactiveValues(
      temp.agregate = NULL
    )
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      file <- paste0(config$path_to_md_dir, '/', id, '.md')
      
      tagList(
        # In this example, the md file is found in the module_examples directory
        # but with a real app, it should be provided by the package which
        # contains the UI for the different steps of the process module.
        # system.file(xxx)
        
        if (file.exists(file))
          includeMarkdown(file)
        else
          p('No Description available'),
        
        
        # Used to show some information about the dataset which is loaded
        # This function must be provided by the package of the process module
        uiOutput(ns('datasetDescription_ui')),
        
        # Insert validation button
        uiOutput(ns('Description_btn_validate_ui'))
      )
    })
    
    output$datasetDescription_ui <- renderUI({
      # Insert your own code to vizualise some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"),
                             "Start",
                             class = btn_success_color)
      toggleWidget(rv$steps.enabled['Description'], widget)
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- Magellan::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- global$VALIDATED
    })
    
    
    
    
    output$Agregation <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by Magellan
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        uiOutput(ns("warningAgregationMethod")),
        
        uiOutput(ns("Agregation_ProteinId_ui")),
        uiOutput(ns('Agregation_useOfShared_ui')),
        uiOutput(ns('Agregation_consider_ui')),
        uiOutput(ns('Agregation_nTopn_ui')),
        uiOutput(ns("Agregation_0perator_ui")),
        
        # Insert validation button
        uiOutput(ns('Agregation_btn_validate_ui')),
        
        # Insert additional code
        uiOutput(ns("ObserverAggregationDone")),
        shinyjs::hidden(downloadButton(ns('downloadAggregationIssues'),
                                       'Download issues', 
                                       class = actionBtnClass)),
        
        hr(),
        uiOutput(ns('Agregation_barplotType_ui')),
        highchartOutput(ns("peptideBarplot"), width="400px"),
        DT::DTOutput(ns("aggregationStats"))
      )
    })
    
    
    output$Agregation_useOfShared_ui <- renderUI({
      mod_popover_for_help_server("modulePopover_useOfShared",
                                  data = list(title="Include shared peptides",
                                              content= HTML(paste0("<ul><li><strong>No:</strong> only protein-specific peptides</li><li><strong>Yes 1:</strong> shared peptides processed as protein specific</li><li><strong>Yes 2</strong>: proportional redistribution of shared peptides</li></ul>")
                                              )
                                  )
      )
      
      widget <- radioButtons(ns("Agregation_useOfShared"), 
                             NULL, 
                             choices = c("No" = "No",
                                         "As protein specific"= "asSpec" ,
                                         "For redistribution" = "forDistrib" ),
                             selected = rv.widgets$Agregation_useOfShared)
      
      tagList(
        mod_popover_for_help_ui(ns("modulePopover_useOfShared")),
        toggleWidget(rv$steps.enabled['Agregation'], widget )
      )

    })

    output$Agregation_ProteinId_ui <- renderUI({
      # req (is.null(GetProteinId(rv$dataIn))
      # 
      # widget <- selectInput(ns("Agregation_ProteinId"), 
      #             "Choose the protein ID",
      #             choices = c("None", colnames(Biobase::fData(rv$dataIn))),
      #             selected = rv.widgets$Agregation_proteinId)
      # toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    
    output$Agregation_consider_ui <- renderUI({
      
      widget <- radioButtons(ns("Agregation_consider"), "Consider",
                             choices=c('all peptides' = "allPeptides", 
                                       "N most abundant" = "onlyN"),
                             selected = rv.widgets$Agregation_consider)
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    
    
    output$Agregation_nTopn_ui <- renderUI({
      req (rv.widgets$Agregation_consider == 'onlyN')
      
      widget <- numericInput(ns("Agregation_nTopn"), "N",
                             value = rv.widgets$Agregation_nTopn, 
                             min = 0, 
                             step = 1, 
                             width = '100px')
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    output$Agregation_operator_ui <- renderUI({
      req(rv.widgets$Agregation_useOfShared)
      
      widget <- radioButtons(ns("Agregation_operator"), "Operator", 
                             choices = if (rv.widgets$Agregation_useOfShared %in% c("No", "asSpec")){
                               choice <- c("Mean"="Mean","Sum"="Sum")
                             } else {choice <- c("Mean"="Mean")}, 
                             selected = rv.widgets$Agregation_operator)
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    
    output$warningAgregationMethod <- renderUI({
      req(rv$dataIn)
      m <- match.qMetadata(qMetadata(rv$dataIn, length(rv$dataIn)), 
                           pattern = "missing", 
                           level = 'peptide')
      
      if (length(which(m)) > 0)
      {
        tags$p(style = "color: red;",
               tags$b('Warning:')," Your dataset contains missing values.
    For better results, you should impute them first")
      }
      
    })
    
    
    
    
    
    
    
    
    output$Agregation_barplotType_ui <- renderUI({
      
      widget <- selectInput(ns("Agregation_barplotType"), "Barplot type",
                             choices=c('All peptides' = "all", 
                                       "Only specific peptides" = "onlySpec",
                                       "Only shared peptides" = "onlyShared"),
                             selected = rv.widgets$Agregation_barplotType,
                            width = '200')
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })

    output$peptideBarplot <- renderHighchart({
      req(adjacencyMatrix(rv$dataIn[[length(rv$dataIn)]]))
     withProgress(message = 'Rendering plot, pleast wait...',detail = '', value = 1, {
       X <- adjacencyMatrix(GetCurrentSE())
       if (is.null(X)){
         X <- makeAdjacencyMatrix(rowData(GetCurrentSE())[,'Protein_group_IDs'])
          rownames(X) <- rownames(GetCurrentSE())
       }
       X <- updateAdjacencyMatrix(X, mode = rv.widgets$Agregation_barplotType)
      GraphPepProt_hc(X)
       })
    })

    
    output$Agregation_nbPeptides_ui <- renderUI({
      req(isTRUE(rv.widgets$Agregation_filterProtAfterAgregation))
      
      widget <- numericInput(ns("Agregation_nbPeptides"), 
                             "Nb of peptides defining a protein", 
                             value = 0, min =0, step=1,
                             width = "250px")
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    
    # output$aggregationStats <- DT::renderDataTable (server=TRUE,{
    #   req(DAPAR::GetMatAdj(rv$current.obj))
    #   req(rv$widgets$aggregation$proteinId != "None")
    #   
    #   res <- getProteinsStats(DAPAR::GetMatAdj(rv$current.obj)$matWithSharedPeptides)
    #   
    #   rv$AggregProtStats$nb <- c(res$nbPeptides,
    #                              res$nbSpecificPeptides,
    #                              res$nbSharedPeptides,
    #                              res$nbProt,
    #                              length(res$protOnlyUniquePep),
    #                              length(res$protOnlySharedPep),
    #                              length(res$protMixPep))
    #   
    #   df <- as.data.frame(rv$AggregProtStats)
    #   names(df) <- c('Description', 'Value')
    #   DT::datatable(df, 
    #                 escape = FALSE,
    #                 rownames= FALSE,
    #                 extensions = c('Scroller'),
    #                 option=list(initComplete = initComplete(),
    #                             dom = 'rt',
    #                             autoWidth = TRUE,
    #                             ordering = F,
    #                             columnDefs = list(list(width='150px', targets= 0),
    #                                               list(width='100px', targets= 1))
    #                 )
    #   )
    # })
    
    
    GetCurrentSE <- reactive({
      rv$dataIn[[length(rv$dataIn)]]
    })
    
    ########################################################
    RunAggregation <- reactive({
      
      # By default, save the whole adjacency matrix
      fcol <- 'Protein_group_IDs'
      X.all <- makeAdjacencyMatrix(rowData(GetCurrentSE())[,fcol])
      rownames(X.all) <- rownames(GetCurrentSE())
      browser()
      rowData(GetCurrentSE())['adjacencyMatrix'] <- NULL
      adjacencyMatrix(GetCurrentSE()) <- X.all
      
      
      
      # withProgress(message = '',detail = '', value = 0, {
      #   incProgress(0.2, detail = 'loading foreach package')
      #   require(foreach)
      #   
      #   incProgress(0.5, detail = 'Aggregation in progress')
      #   
         ll.agg <- NULL
         
         if(rv.widgets$Agregation_useOfShared == "forDistrib"){
           if (rv.widgets$Agregation_consider == 'allPeptides'){
             # aggIterative agregation method
             ll.agg <- aggregateFeatures4Prostar(object = rv$dataIn,
                                                 i = length(rv$dataIn), 
                                                 name = 'aggregated',
                                                 fcol = "adjacencyMatrix", 
                                                 fun = aggIterative,
                                                 conditions = colData(object)$Condition, 
                                                 init.method = 'rowSums',
                                                 iter.method = 'Mean')
           } else if (rv.widgets$Agregation_consider == 'onlyN'){
             ll.agg <- aggregateFeatures4Prostar(object = rv$dataIn,
                                                 i = length(rv$dataIn), 
                                                 name = 'aggregated',
                                                 fcol = "adjacencyMatrix", 
                                                 fun = aggIterative,
                                                 conditions = colData(object)$Condition, 
                                                 init.method = 'rowSums',
                                                 iter.method = 'onlyN',
                                                 n = 10)
             }
         } else if (rv.widgets$Agregation_useOfShared == "asSpec") {
           if (rv.widgets$Agregation_consider == 'allPeptides') {
             X.all <- X.all
           } else if (rv.widgets$Agregation_consider == 'onlyN') {
             rowData(GetCurrentSE())[['adjacencyMatrix']] <- NULL
             X.n <- updateAdjacencyMatrix(X,
                                          mode = 'topn',
                                          qData = assay(GetCurrentSE()),
                                          fun = 'rowMeans',
                                          n = rv.widgets$Agregation_topn)
             adjacencyMatrix(GetCurrentSE()) <- X.n
             }
           ll.agg <- aggregateFeatures4Prostar(object = rv$dataIn,
                                               i = length(rv$dataIn), 
                                               name = 'aggregated',
                                               fcol = "adjacencyMatrix", 
                                               fun = rv.widgets$Agregation_operator)
           
           
         } else if (rv.widgets$Agregation_useOfShared == "no") {
           
           # Extract only specific peptides for adjacency matrix
           rowData(GetCurrentSE())[['adjacencyMatrix']] <- NULL
           X.spec <- updateAdjacencyMatrix(adjacencyMatrix(GetCurrentSE()),
                                        mode = 'onlySpec',
                                        qData = assay(GetCurrentSE())
                                        )
           adjacencyMatrix(GetCurrentSE()) <- X.spec
           
           if (rv.widgets$Agregation_consider == 'onlyN') {
             # Extract only n best peptides for adjacency matrix
             rowData(GetCurrentSE())[['adjacencyMatrix']] <- NULL
             X.best.n <- updateAdjacencyMatrix(adjacencyMatrix(GetCurrentSE()),
                                          mode = 'topn',
                                          qData = assay(GetCurrentSE()),
                                          fun = 'rowMeans',
                                          n = rv.widgets$Agregation_topn)
             adjacencyMatrix(GetCurrentSE()) <- X.best.n
           }
           ll.agg <- aggregateFeatures4Prostar(object = rv$dataIn,
                                               i = length(rv$dataIn), 
                                               name = 'aggregated',
                                               fcol = "adjacencyMatrix", 
                                               fun = rv.widgets$Agregation_operator)
         }
      return(ll.agg)
      
    })

    
    
    output$displayNbPeptides <- renderUI({
      req(rv$widgets$aggregation$filterProtAfterAgregation)
      
      if (rv$widgets$aggregation$filterProtAfterAgregation) {
        numericInput("nbPeptides", "Nb of peptides defining a protein", 
                     value = 0, min =0, step=1,
                     width = "250px")
      }
    })
    
    output$Agregation_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Agregation_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      toggleWidget(rv$steps.enabled['Agregation'], widget)
    })
    
    
    observeEvent(input$Agregation_btn_validate, {
      print('Perform')
     # rv$dataIn <- RunAggregation()
      
      
      
      rv.custom$temp.agregate <- RunAggregation()
      
      if(is.null(rv.custom$temp.agregate$issues)){
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
        rv$steps.status['Agregation'] <- global$VALIDATED
      } else {
        print("error")
        shinyjs::toggle('downloadAggregationIssues', 
                        condition = !rvModProcess$moduleAggregationDone[1] && length(rv.custom$temp.agregate$issues) > 0
        )
      }
      

    })
    

    
    output$downloadAggregationIssues <- downloadHandler(
      
      filename = 'aggregation_issues.txt',
      content = function(file) {
        
        tmp.peptides <- lapply(rv.custom$temp.agregate$issues, function(x)paste0(x, collapse=","))
        df <- data.frame(Proteins=names(rv.custom$temp.agregate$issues), 
                         Peptides = as.data.frame(do.call(rbind, tmp.peptides))
        )
        colnames(df) <- c('Proteins', 'Peptides')
        write.table(df, file = file, row.names = FALSE, quote=FALSE, sep="\t")
      }
    )
    
    
    
    
    #--------------------------------------------------------
    output$Save <- renderUI({
      wellPanel(
        uiOutput(ns('AddMetadata_btn_validate_ui'))
      )
    })
    
    output$AddMetadata_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Save_btn_validate"),
                             "Perform",
                             class = btn_success_color)
      toggleWidget(rv$steps.enabled['Save'], widget)
    })
    
    
    observeEvent(input$Save_btn_validate, {
      
      
      
      #dataOut$trigger <- Magellan::Timestamp()
      #dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- global$VALIDATED
    })
    
    
    #----------------------------------------------------
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

