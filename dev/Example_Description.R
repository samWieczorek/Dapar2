#Timeline_R6.R
Example_Description = R6::R6Class(
  "Example_Description",
  inherit = Process,
  private = list(
    .config = list(name = 'Description',
                   steps = c('Description'),
                   mandatory = c(T)
    )
  ),

  public = list(

    Global_server = function(input, output){},


    Description_server = function(input, output){
      observeEvent(input$btn_validate_Description, ignoreInit = T, ignoreNULL=T, {
        cat(paste0(class(self)[1], "::observeEvent(input$btn_validate_Description from - ", self$id, '\n'))
        private$InitializeDataIn()
        self$ValidateCurrentPos()
      })

      output$datasetDescription <- renderUI({
        tagList(
          p(paste0('Dataset description: ', paste0(names(self$rv$temp.dataIn), collapse=", ")))
        )
      })
    },


    Description_ui = function(){
      wellPanel(
        tagList(
        includeMarkdown( system.file("app/md", paste0(self$config$name, ".md"), package="Magellan")),
        uiOutput(self$ns('datasetDescription')),
        actionButton(self$ns('btn_validate_Description'),
                     paste0('Start ', self$config$name),
                     class = 'btn-success')
      )
      )
    }
  )
)
