
source("./chunks_in_variables.R") #names(.GlobalEnv)


write_ui <- function(name, file) {
  new_ui <- gsub('#name#', name, ui)
  cat(new_ui, file=file, append = TRUE)
}


write_start_server <- function(name, file) {
  new_start_server <- gsub('#name#' ,name , start_server)
  cat(new_start_server, file=file)
}


write_widgets_default_values <- function(name, widgets.names, widgets.defaultVal, file){
  cat("# Define default selected values for widgets\n", file=file, append = TRUE)
  cat("widgets.default.values <- list(\n", file=file, append = TRUE)
  if (length(widgets.names) > 0)
      for (i in 1:length(widgets.names)){
        cat(paste0(widgets.names[i], " = ", widgets.defaultVal[i]), file=file, append = TRUE)
        if (i < length(widgets.names))
          cat(",", file=file, append = TRUE)
        cat("\n", file=file, append = TRUE)
  }

  cat(")\n\n", file=file, append = TRUE)
}


write_moduleServerFunc <- function(file){
  cat(moduleServerFunc_template, file=file, append = TRUE)
}

write_config <- function(name, steps, mandatory, file) {
  new_config <- gsub('#name#' ,name , config_template)
  new_config <- gsub('#steps#' ,steps , new_config)
  new_config <- gsub('#mandatory#' ,mandatory , new_config)
  cat(new_config, file=file, append = TRUE)
}


write_rv_widgets_def <- function(widgets.names, file){
  if (length(widgets.names)==0) return(NULL)
  cat("rv.widgets <- reactiveValues(\n", file=file, append = TRUE)

  for (i in 1:length(widgets.names)){
    cat(paste0(widgets.names[i], " =  widgets.default.values$", widgets.names[i]), file=file, append = TRUE)
    if (i < length(widgets.names))
      cat(",", file=file, append = TRUE)
    cat("\n", file=file, append = TRUE)
  }

  cat(")\n\n", file=file, append = TRUE)
}



write_initModule_func <- function(file){
  cat(chunk_init_module_template, file=file, append = TRUE)
}

Add_observeEvent_stepsEnabled <- function(file){
  cat(chunk_observeEvent_steps_enabled, file=file, append = TRUE)
}

Add_observeEvent_remoteReset <- function(file){
  cat(chunk_observeEvent_remoteReset, file=file, append = TRUE)
}


write_output_description_template <- function(file){
  cat(output_Description_template, file=file, append = TRUE)

}


write_header_comment <- function(step.name, file){
  txt <- gsub('#step.name#', step.name, write_comment_for_code)
  cat(txt, file=file, append = TRUE)
}


Add_ObserveEvent_rv_comment <- function(file){
  cat(observer_for_widget_comment, file=file,append=TRUE)
}

write_observer_for_widgets <- function(widgets.names, file){
  if (length(widgets.names)>0)
  for (i in 1:length(widgets.names))
    cat(gsub('#name#' , widgets.names[i] , observer_for_widget_template),
        file=file, append = TRUE)
}


write_code_for_renderUI_widgets <- function(widgets.renderUI.names, step.name, file){
  txt <- gsub('#step#' , step.name , widget_renderUI_template)

  for (i in widgets.renderUI.names)
    cat(gsub('#name#' , i , txt), file = file, append = TRUE)
}



write_header_for_global_step_renderUI <- function(step.name, file){
  cat(gsub('#step.name#' , step.name , code_for_global_step_renderUI_header),
      file = file,
      append = TRUE)
}

write_code_for_call_renderUI_widgets <- function(widgets.renderUI.names, file){

  for (i in 1:length(widgets.renderUI.names)){
    cat(gsub('#widget.name#' , widgets.renderUI.names[i] , code_for_global_step_renderUI_call_widget_renderUI),
        file = file, append = TRUE)
    if (i < length(widgets.renderUI.names))
      cat(',', file = file, append = TRUE)

  }
}

write_additional_comma <- function(file){
  cat(',', file = file, append = TRUE)
}

write_code_for_call_direct_widgets <- function(widgets.direct.names, step.name, file){
  txt <- gsub('#step.name#' , step.name , code_for_global_step_renderUI_call_widget)
  for (i in 1:length(widgets.direct.names)){
    cat(gsub('#widget.name#', widgets.direct.names[i], txt),
        file = file, append = TRUE)
    if (i < length(widgets.direct.names))
      cat(',', file = file, append = TRUE)

  }
}

write_code_for_validation_btn <- function(step.name, file){
  txt <- gsub('#step.name#' , step.name , code_for_global_step_renderUI_call_validate_btn)
  cat(txt, file =file, append = TRUE)
}


write_code_for_observe_validation_btn_generic_step <- function(step.name, file){

  txt <- gsub('#step.name#' , step.name , chunk_for_observe_validation_btn_generic_step)
  cat(txt, file =file, append = TRUE)
}

write_code_for_observe_validation_btn_last_step <- function(step.name, file){

  txt <- gsub('#step.name#' , step.name , chunk_for_observe_validation_btn_last_step)
  cat(txt, file =file, append = TRUE)
}



write_code_to_end_global_step_renderUI <- function(file){

  cat(code_for_global_step_renderUI_end, file = file, append = TRUE)
}


add_code_for_validation_btn_generic_step <- function(step.name, file){
  txt <- gsub('#step.name#' , step.name , chunk_for_validation_btn_generic_step)
  cat(txt, file =file, append = TRUE)
}


add_code_for_validation_btn_last_step <- function(step.name, file){
  txt <- gsub('#step.name#' , step.name , chunk_for_validation_btn_last_step)
  cat(txt, file =file, append = TRUE)
}


write_code_for_return_server <- function(file){
  cat(chunk_for_return_server, file = file, append = TRUE)
}

write_code_for_end_module <- function(file){
  cat(chunk_for_end_module, file = file, append = TRUE)
}
