library(Magellan)
library(QFeatures)
library(DaparToolshed)

#source(file.path('.', 'global.R'), local=TRUE)$value
#source(file.path('.', 'mod_Agregation.R'), local=TRUE)$value

options(shiny.fullstacktrace = FALSE,
        shiny.maxRequestSize=30*1024^2) 

run_workflow('Agregation', verbose=T)