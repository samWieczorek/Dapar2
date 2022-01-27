
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
    steps = c('Description', 'Agregation', 'Add metadata', 'Save'),
    
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, TRUE, FALSE, TRUE),
    
    path_to_md_dir =  system.file('md/', package='DaparToolshed')
  )
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    Agregation_ProteinId = NULL,
    Agregation_includeShared = NULL,
    Agregation_Consider = NULL,
    Agregation_nTopn = NULL,
    Agregation_Operator = NULL,
    nbPeptides = 0,
    Agregation_filterProtAfterAgregation = FALSE
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
        uiOutput(ns('Agregation_includeShared_ui')),
        uiOutput(ns('Agregation_Consider_ui')),
        uiOutput(ns('Agregation_nTopn_ui')),
        uiOutput(ns("Agregation_0perator_ui")),
        
        # Insert validation button
        uiOutput(ns('Agregation_btn_validate')),
        
        # Insert additional code
        uiOutput(ns("ObserverAggregationDone")),
        shinyjs::hidden(downloadButton(ns('downloadAggregationIssues'),
                                       'Download issues', 
                                       class = actionBtnClass)),
        
        hr(),
        div(
          div( style="display:inline-block; vertical-align: top;",
               uiOutput(ns("specificPeptideBarplot"))),
          div( style="display:inline-block; vertical-align: top; padding-right: 20px;",       
               uiOutput(ns("allPeptideBarplot"))),
          div( style="display:inline-block; vertical-align: top;",
               tagList(
                 DT::DTOutput(ns("aggregationStats"))
               )
          )
        )
        
      )
    })
    
    
    output$Agregation_includeShared_ui <- renderUI({
      mod_popover_for_help_server("modulePopover_includeShared",
                                  data = list(title="Include shared peptides",
                                              content= HTML(paste0("<ul><li><strong>No:</strong> only protein-specific peptides</li><li><strong>Yes 1:</strong> shared peptides processed as protein specific</li><li><strong>Yes 2</strong>: proportional redistribution of shared peptides</li></ul>")
                                              )
                                  )
      )
      
      widget <- radioButtons(ns("Agregation_includeShared"), 
                             NULL, 
                             choices = c("No" = "No",
                                         "Yes (as protein specific)"= "Yes1" ,
                                         "Yes (redistribution)" = "Yes2" ),
                             selected = rv.widgets$Agregation_includeShared)
      
      tagList(
        mod_popover_for_help_ui(ns("modulePopover_includeShared")),
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
    
    
    
    output$Agregation_Consider_ui <- renderUI({
      
      widget <- radioButtons("Agregation_Consider", "Consider",
                             choices=c('all peptides' = "allPeptides", 
                                       "N most abundant" = "onlyN"),
                             selected = rv.widgets$Agregation_Consider)
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    
    
    output$Agregation_nTopn_ui <- renderUI({
      browser()
      req (rv.widgets$Agregation_Consider == 'onlyN')
      
      widget <- numericInput(ns("Agregation_nTopn"), "N",
                             value = rv.widgets$Agregation_nTopn, 
                             min = 0, 
                             step = 1, 
                             width = '100px')
      
      toggleWidget(rv$steps.enabled['Agregation'], widget )
    })
    
    
    output$Agregation_Operator_ui <- renderUI({
      req(rv.widgets$Agregation_includeShared)
      
      # choice <- NULL
      # if (rv.widgets$Agregation_includeShared %in% c("No", "Yes1")){
      #   choice <- c("Mean"="Mean","Sum"="Sum")
      # } else {choice <- c("Mean"="Mean")}
      # choice
      
      widget <- radioButtons(ns("Agregation_operator"), "Operator", 
                             choices = if (rv.widgets$Agregation_includeShared %in% c("No", "Yes1")){
                               choice <- c("Mean"="Mean","Sum"="Sum")
                             } else {choice <- c("Mean"="Mean")}, 
                             selected = rv.widgets$Agregation_Operator)
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
    
    observeEvent(rv.widgets$Agregation_includeShared, {
      if (rv.widgets$Agregation_includeShared == 'Yes2'){
        ch <- c("Mean" = "Mean")  
      } else {
        ch <- c("Sum" = 'Sum', "Mean"="Mean")
      }
    })
    
    
    
    
    output$specificPeptideBarplot <- renderUI({
      req(GetAdjMat())
      withProgress(message = 'Rendering plot, pleast wait...',detail = '', value = 1, {
        tagList(
          h4("Only specific peptides"),
          plotOutput(ns("aggregationPlotUnique"), width="400px")
        )
      })
    })
    
    output$allPeptideBarplot <- renderUI({
      req(GetAdjMat())
      withProgress(message = 'Rendering plot, pleast wait...',detail = '', value = 1, {
        tagList(
          h4("All (specific & shared) peptides"),
          plotOutput(ns("aggregationPlotShared"), width="400px")
        )
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
    
    
    GetAdjMat <- reactive({
      #browser()
      req(rv$dataIn)
      rowData(rv$dataIn[[length(rv$dataIn)]])$adjacencyMatrix
    })
    
    ########################################################
    RunAggregation <- reactive({
      req(GetAdjMat())
      rv.widgets$Agregation_includeShared
      rv.widgets$Agregation_operator
      rv.widgets$Agregation_consider
      rv.widgets$Agregation_topN
      
      withProgress(message = '',detail = '', value = 0, {
        incProgress(0.2, detail = 'loading foreach package')
        
        
        require(foreach)
        incProgress(0.5, detail = 'Aggregation in progress')
        
        ll.agg <- NULL
        if(rv.widgets$Agregation_includeSharedPeptides %in% c("Yes2", "Yes1")){
          X <- submatadj(GetAdjMat(), onlyShared = TRUE)
          if (rv.widgets$Agregation_includeSharedPeptides == 'Yes1'){
            if (rv.widgets$Agregation_consider == 'allPeptides') {
              ll.agg <- do.call(paste0('aggregate', rv.widgets$Agregation_operator),
                                list( obj.pep = rv$dataIn, 
                                      X = X))
            } else {
              ll.agg <- aggregateTopn(rv$dataIn, 
                                      X,
                                      rv.widgets$Agregation_operator, 
                                      n = as.numeric(rv.widgets$Agregation_topN))
            }
          } else {
            if (rv.widgets$Agregation_consider == 'allPeptides') {
              ll.agg <- aggregateIterParallel(obj.pep = rv$dataIn, 
                                              X = X,
                                              init.method = 'Sum', 
                                              method = 'Mean')
            } else {
              ll.agg <- aggregateIterParallel(rv$dataIn, 
                                              X, 
                                              init.method = 'Sum', 
                                              method = 'onlyN', 
                                              n = rv.widgets$Agregation_topN)
            }
          }
        } else {
          X <- submatadj(GetAdjMat(), onlySpec = TRUE)
          if (rv.widgets$Agregation_consider == 'allPeptides') {
            ll.agg <- do.call(paste0('aggregate', rv.widgets$Agregation_operator),
                              list(obj.pep = rv$dataIn,
                                   X = X))
          } else {
            ll.agg <- aggregateTopn(rv$dataIn, 
                                    X = X, 
                                    rv.widgets$Agregation_operator,
                                    n = as.numeric(rv.widgets$Agregation_topN)
            )
          }
        }
      } )
      
      
      
      return(ll.agg)
      
    })
    
    
    
    
    ###------------ Perform aggregation--------------------
    observeEvent(input$Agregation_btn_validate,{
      
      rv.custom$temp.agregate <- RunAggregation()
      
      if(is.null(rv.custom$temp.agregate$issues)){
        dataOut$trigger <- Magellan::Timestamp()
        dataOut$value <- rv$dataIn
        rv$steps.status['Step1'] <- global$VALIDATED
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
    
    
    #-------------------------------------------------------
    
    output$Addmetadata <- renderUI({ })
    
    
    output$Save <- renderUI({ })
    
    
    
    # Insert necessary code which is hosted by Magellan
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

