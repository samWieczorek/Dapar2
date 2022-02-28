

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
#' 
#' @examples 
#' create_ft_example()
#' create_ft_example(TRUE)
#' 
create_ft_example <- function(with.na = FALSE){
  
  filename <- if (with.na) "ft-data-na.txt" else "ft-data.txt"
  data.file <- system.file("extdata", filename, package="DaparToolshed")
  data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
  sample.file <- system.file("extdata", "ft-samples.txt", package="DaparToolshed")
  sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)

  tmp.qf <- createQFeatures(data = data,
                            sample = sample,
                            indQData = 2:7,
                            keyId = 'ID',
                            analysis = 'test',
                            indQMetadata = 9:14,
                            typeDataset = 'peptide',
                            parentProtId = 'Proteins',
                            force.na = TRUE,
                            software = 'maxquant')

  
  if(with.na){
    ft_na <- tmp.qf
    save(ft_na, file= paste0('data/ft_na.rda'))
  } else {
    ft <- tmp.qf
    save(ft, file= paste0('data/ft.rda'))
  }
  
  tmp.qf
  }