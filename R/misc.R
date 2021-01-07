#' @title
#' xxx
#' 
#' @description
#' Builds a vector of `#conditions` colors. 
#' 
#' @param conditions xxxx
#' 
#' @param base_palette xxx
#' 
#' @export
#' 
Example_Palette <- function(conditions, base_palette){
  
  examplePalette <- NULL
  nbConds <- length(unique(conditions))
  for (i in 1:nbConds){
    examplePalette[ which(conditions == unique(conditions)[i])] <- base_palette[i]
  }
  
  return(examplePalette)
}



#' @title
#' xxx
#' 
#' @description
#' Builds a vector of colors (one per condition) to be used as a basis palette
#' 
#' @param palette.name The name of the palette from the package `RColorBrewer`
#' 
#' @param conditions xxxx
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#'
Base_Palette <- function(palette.name = 'Dark2', conditions){

nbConds <- length(unique(conditions))
nbColors <- max(3, nbConds)
basePalette <- RColorBrewer::brewer.pal(nbColors, palette.name)[1:nbConds]
return(basePalette)
}