library(shiny)

options(shiny.fullstacktrace = TRUE)
options(shiny.trace=FALSE)

source("./process_template_funcs.R")
source("./mod_AddStep.R")




# download button


ui <- fluidPage(
  tagList(
    textInput('pipeline_name','Pipeline name', value = 'pipe1'),
  textInput('process_name','process name', value = 'proc1'),
  selectInput("members", "Number of steps", choices = 2:10, width='100px'),
  uiOutput('show_step_forms'),

  actionButton('validate', 'Generate process skeleton')
  )
)



server <- function(input, output, session){

  rv <- reactiveValues(
    step_form = list(),
    step_value = NULL
  )

  output$show_step_forms <- renderUI({
    members <- as.integer(input$members) # default 2
    req(input$members > 0)
    rv$step_value <- lapply(1:as.numeric(input$members), function(i) {
      mod_AddStep_server(paste0('step_', i),
                         n = reactive({i}),
                         nTotal = reactive({as.numeric(input$members)}))
    })
    lapply(1:members, function(i) {
      mod_AddStep_ui(paste0('step_', i))
    })
  })


  observeEvent(input$validate, {
    call_createModule()
    })

  call_createModule <- function(){
   # browser()

    pipeline <- paste0('mod_', input$pipeline_name)
    process <- input$process_name
    #----------------------------------------------------------------#


    file <- paste0("./test/temp_mod_pipe_", process, ".R")

    cat(NULL, file=file)


    #----------------------------------------------------------------#

    ll.steps <- lapply(1:length(rv$step_value), function(x){rv$step_value[[x]]$name()})
    vec.steps <- unlist(ll.steps) #separe chaque etape
    steps <- paste0('\'',vec.steps,'\'') #chaque mot encadre
    steps <- paste0(steps,collapse = ",") #chaine char sep par ,
    steps <- paste0('c(', steps, ')') #chaine char like vector

    #----------------------------------------------------------------#

    ll.isMandatory <-  lapply(1:length(rv$step_value), function(x){rv$step_value[[x]]$isMandatory()})
    mandatory <- unlist(ll.isMandatory) #separe chaque etape
    mandatory <- paste0(mandatory,collapse = ",") #chaine char sep par ,
    mandatory <- paste0('c(',mandatory,')') #chaine char like vector

    #----------------------------------------------------------------#



    write_ui(name=process,
             file=file)

    write_start_server(name=process,
                       file=file)

    all.widgets.name <- unlist(lapply(1:length(rv$step_value), function(x){rv$step_value[[x]]$widgets.name()}))
    all.widgets.defaultVal <- unlist(lapply(1:length(rv$step_value), function(x){rv$step_value[[x]]$widgets.defaultVal()}))
    write_widgets_default_values(name = process,
                                 widgets.name = all.widgets.name,
                                 widgets.defaultVal = all.widgets.defaultVal,
                                 file = file)

    write_moduleServerFunc(file = file)


    write_config(name = process,
                 steps = steps,
                 mandatory = mandatory,
                 file = file)


    write_rv_widgets_def(widgets.name = all.widgets.name,
                         file = file)


    write_initModule_func(file = file)
    Add_observeEvent_stepsEnabled(file)
    Add_observeEvent_remoteReset(file)

    write_output_description_template(file = file)

    # Add code for each step except the first one
    for (i in 2:length(vec.steps)){
      step.widgets.name <- unlist(rv$step_value[[i]]$widgets.name())
      step.widgets.defaultVal <- unlist(rv$step_value[[i]]$widgets.defaultVal())
      step.widgets.renderUI <- unlist(rv$step_value[[i]]$widgets.renderUI())
      step.widgets.type <- unlist(rv$step_value[[i]]$widgets.type())

      write_header_comment(step.name = vec.steps[i], file = file)
     #browser()
     if (length(step.widgets.name) > 0){
       Add_ObserveEvent_rv_comment(file)
       write_observer_for_widgets(step.widgets.name, file)
     }

      ind.renderUI.widgets <- which(step.widgets.renderUI == TRUE)
      if (length(ind.renderUI.widgets) > 0){
        write_code_for_renderUI_widgets(widgets.renderUI.names = step.widgets.name[ind.renderUI.widgets],
                                        step.name = vec.steps[i],
                                        file = file)

      }

      write_header_for_global_step_renderUI(step.name = vec.steps[i], file = file)

      if (length(ind.renderUI.widgets) > 0){
        write_code_for_call_renderUI_widgets(widgets.renderUI.names = step.widgets.name[ind.renderUI.widgets],
                                             file = file)
      write_additional_comma(file)
      }

      ind.direct.widgets <- which(step.widgets.renderUI == FALSE)
      if (length(ind.direct.widgets) > 0){
        write_code_for_call_direct_widgets(widgets.direct.names = step.widgets.name[ind.direct.widgets],
                                           step.name =  vec.steps[i],
                                           file = file)
        write_additional_comma(file)
      }

      write_code_for_validation_btn(step.name = vec.steps[i], file = file)
      write_code_to_end_global_step_renderUI(file = file)

      if (i == length(vec.steps))
        write_code_for_observe_validation_btn_last_step(step.name = vec.steps[i], file = file)
      else
        write_code_for_observe_validation_btn_generic_step(step.name = vec.steps[i], file = file)



    }
    write_code_for_return_server(file = file)
    write_code_for_end_module(file = file)

  }

}


shinyApp(ui, server)
