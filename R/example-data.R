

#' @title Creates an object of class [QFeatures] as example
#' 
#' @description 
#' xxxxx
#' 
#' @return An example object of class [QFeatures]
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' @rdname create-example
create_ft_example <- function(){
  data.file <- system.file("extdata", "ft-data.txt", package="DaparToolshed")
  data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
  sample.file <- system.file("extdata", "ft-samples.txt", package="DaparToolshed")
  sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)

  ft <- createQFeatures(data = data,
                        sample = sample,
                        indQData = 2:7,
                        keyId = 'ID',
                        analysis = 'test',
                        indQMetadata = 9:14,
                        typeDataset = 'peptide',
                        parentProtId = 'Proteins',
                        force.na = TRUE,
                        software = 'maxquant')

   save(ft, file='data/ft.rda')
  ft
  }