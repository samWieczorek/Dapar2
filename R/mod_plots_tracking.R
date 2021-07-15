#' @title   mod_plots_boxplots_ui and mod_plots_boxplots_server
#'
#' @description  A shiny Module.
#'
#' @param id shiny id
#'
#' @param input internal
#'
#' @param output internal
#'
#' @param session internal
#'
#' @rdname mod_plots_tracking
#'
#' @keywords internal
#'
#' @export
#'
#' @importFrom shiny NS tagList
#'
mod_plots_tracking_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    # shinyjs::hidden(actionButton(ns('rst_btn'), 'Reset')),
    #uiOutput(ns('selectType_UI')),
    selectInput(ns("typeSelect"), "Type of selection",
                choices=c("None" = "None",
                          "Protein list" = "ProteinList",
                          "Random" = "Random",
                          "Specific column" = "Column"),
                width=('130px')),
    shinyjs::hidden(
      selectInput(ns("listSelect"),
                  "Protein for normalization",
                choices = c('None'),
                  multiple = TRUE,
                  width='400px'
      )
      ),
    shinyjs::hidden(
      textInput(ns("randSelect"),
                "Random",
                width = ('120px'))
      ),
    shinyjs::hidden(
      selectInput(ns("colSelect"),
                  "Column",
                  choices = c('')
                )
    )
  )
}

#' plots_tracking Server Function
#'
#' @param obj Object SummarizedExperiment
#'
#' @param metadata Metadata of Features containing the SummarizedExperiment
#'
#' @rdname mod_plots_tracking
#'
#' @export
#'
#' @keywords internal
#'
#' @import shinyjs
#'
mod_plots_tracking_server <- function(id,
                                      obj,
                                      keyId,
                                      params = NULL,
                                      reset = FALSE,
                                      slave = reactive({FALSE})){

  moduleServer(id, function(input, output, session){
    ns <- session$ns

    rv.track <- reactiveValues(
      res = list(typeSelect = "None",
                 listSelect = NULL,
                 randSelect = '',
                 colSelect = NULL,
                 list.indices = NULL,
                 rand.indices = NULL,
                 col.indices = NULL),
      sync = FALSE
    )


    observe({
      req(obj())
      if (class(obj()) != "SummarizedExperiment") { return(NULL) }

      rv.track$res <- list(typeSelect = if (is.null(params()$typeSelect)) 'None' else params()$typeSelect,
                           listSelect = params()$listSelect,
                           randSelect = params()$randSelect,
                           colSelect = params()$colSelect,
                           list.indices = NULL,
                           rand.indices = '',
                           col.indices = NULL)


      updateSelectInput(session, 'typeSelect',
                        selected = rv.track$res$typeSelect)
      updateSelectInput(session, 'listSelect',
                        choices = SummarizedExperiment::rowData(obj())[[keyId()]],
                        selected = rv.track$res$listSelect)
      updateTextInput(session, 'randSelect',
                      value = rv.track$res$randSelect)
      updateSelectInput(session, 'colSelect',
                        selected = rv.track$res$colSelect,
                        colnames(SummarizedExperiment::rowData(obj()))
      )

      shinyjs::toggle(id = 'listSelect', condition = input$typeSelect == 'ProteinList')
      shinyjs::toggle(id = 'randSelect', condition = input$typeSelect == 'Random')
      shinyjs::toggle(id = 'colSelect', condition = input$typeSelect == 'Column')

    })


    observeEvent(req(input$typeSelect),{

      rv.track$res <- list(typeSelect = if (is.null(input$typeSelect)) 'None' else input$typeSelect,
                           listSelect = NULL,
                           randSelect = '',
                           colSelect = NULL,
                           list.indices = NULL,
                           rand.indices = '',
                           col.indices = NULL)
    })



    observeEvent(slave(),{
      if(is.null(slave()))
        rv.track$sync <- FALSE
      else
        rv.track$sync <- slave()
    })





    # output$selectType_UI <- renderUI({
    #   selectInput(ns("typeSelect"), "Type of selection",
    #               choices=c("None" = "None",
    #                         "Protein list" = "ProteinList",
    #                         "Random" = "Random",
    #                         "Specific column" = "Column"),
    #               selected = params()$typeSelect,
    #               width=('130px'))
    # })


    # output$listSelect_UI <- renderUI({
    #   selectInput(ns("listSelect"),
    #               "Protein for normalization",
    #               choices = SummarizedExperiment::rowData(obj())[[keyId()]],
    #               selected = params()$listSelect,
    #               multiple = TRUE,
    #               width='400px'
    #   )
    # })
    #
    #
    # output$randomSelect_UI <- renderUI({
    #   textInput(ns("randSelect"),
    #             "Random",
    #             value = params()$randSelect,
    #             width = ('120px'))
    # })
    #
    # output$columnSelect_UI <- renderUI({
    #   selectInput(ns("colSelect"),
    #               "Column",
    #               choices = colnames(SummarizedExperiment::rowData(obj())),
    #               selected = params()$colSelect
    #               )
    # })




    observe({

      reset()
      if (reset() > 0) {
        updateSelectInput(session, "typeSelect", selected='None')
        updateSelectInput(session, "listSelect", NULL)
        updateSelectInput(session, "randSelect", selected='')
        updateSelectInput(session, "colSelect", selected=NULL)
        rv.track$res <-list(typeSelect = "None",
                            listSelect = NULL,
                            randSelect = '',
                            colSelect = NULL,
                            list.indices = NULL,
                            rand.indices = NULL,
                            col.indices = NULL)
      }
    })






    observeEvent(rv.track$sync, ignoreNULL = TRUE,{

      if (rv.track$sync == FALSE) {
        rv.track$res <- list(typeSelect = 'None',
                             listSelect = NULL,
                             randSelect = '',
                             colSelect = NULL,
                             list.indices = NULL,
                             rand.indices = '',
                             col.indices = NULL)
        shinyjs::show("typeSelect")
      } else {
        updateSelectInput(session, "typeSelect", selected='None')
        updateSelectInput(session, "listSelect", NULL)
        updateSelectInput(session, "randSelect", selected='')
        updateSelectInput(session, "colSelect", selected=NULL)

        shinyjs::hide("typeSelect")
        shinyjs::hide("listSelect_UI")
        shinyjs::hide("randomSelect_UI")
        shinyjs::hide("columnSelect_UI")
      }

    })









    observeEvent(input$listSelect, ignoreNULL = FALSE,{
      rv.track$res$listSelect <- input$listSelect
      updateSelectInput(session, "randSelect", selected='')
      updateSelectInput(session, "colSelect", selected=NULL)

      if(is.null(rv.track$res$listSelect))
        rv.track$res$list.indices <- NULL
      else
        rv.track$res$list.indices <-  match(rv.track$res$listSelect, SummarizedExperiment::rowData(obj())[[keyId()]])
    })





    observeEvent(input$randSelect,ignoreNULL = FALSE,{
      rv.track$res$randSelect <- input$randSelect

      updateSelectInput(session, "listSelect", NULL)
      updateSelectInput(session, "colSelect", selected=NULL)

      if (is.null(rv.track$res$randSelect) || rv.track$res$randSelect==''
          || (as.numeric(rv.track$res$randSelect) < 0))
      {
        rv.track$res$rand.indices <- NULL
      } else {
        rv.track$res$rand.indices <- sample(1:nrow(obj()), as.numeric(rv.track$res$randSelect), replace=FALSE)
      }
    })

    observeEvent(input$colSelect, ignoreNULL = FALSE,{
      req(obj())
      rv.track$res$colSelect <- input$colSelect

      updateSelectInput(session, "listSelect", NULL)
      updateSelectInput(session, "randSelect", selected='')

      if (is.null(rv.track$res$colSelect))
        rv.track$res$col.indices <- NULL
      else
        rv.track$res$col.indices <- which(SummarizedExperiment::rowData(obj())[[rv.track$res$colSelect]] == 1)

     })

    return(reactive({rv.track$res}))


  })

}

## To be copied in the UI
# mod_plots_tracking_ui("plots_tracking_ui_1")

## To be copied in the server
# callModule(mod_plots_tracking_server, "plots_tracking_ui_1")

