#' #' @title
#' #' xxx
#' #' 
#' #' @description
#' #' Builds a vector of `#conditions` colors. 
#' #' 
#' #' @param conditions xxxx
#' #' 
#' #' @param base_palette xxx
#' #' 
#' #' @export
#' #' 
#' #' @return NA
#' #' 
#' Example_Palette <- function(conditions, base_palette){
#'   
#'   examplePalette <- NULL
#'   nbConds <- length(unique(conditions))
#'   for (i in 1:nbConds){
#'     examplePalette[ which(conditions == unique(conditions)[i])] <- base_palette[i]
#'   }
#'   
#'   return(examplePalette)
#' }
#' 
#' 
#' 
#' #' @title Extends a base-palette of the package RColorBrewer to n colors.
#' #' 
#' #' @description The colors in the returned palette are always in the same order
#' #' 
#' #' @param n The number of desired colors in the palette
#' #' 
#' #' @param base The name of the palette of the package RColorBrewer from which 
#' #' the extended palette is built. Default value is 'Set1'.
#' #' 
#' #' @return A vector composed of n color code.
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' ExtendPalette(12)
#' #' nPalette <- 10
#' #' par(mfrow=c(nPalette,1))
#' #' par(mar=c(0.5, 4.5, 0.5, 0.5))
#' #' for (i in 1:nPalette){
#' #'   pal <- ExtendPalette(n=i, base = 'Dark2')
#' #'   barplot(1:length(pal), col=pal)
#' #'   print(pal)
#' #' }
#' #' 
#' #' @export
#' #' 
#' #' @importFrom RColorBrewer brewer.pal
#' #' @importFrom grDevices colorRampPalette
#' #' 
#' #' @return NA
#' #' 
#' ExtendPalette <- function(n = NULL, base = "Set1"){
#' 
#'   palette <- NULL
#'   if(is.null(base))
#'     base <- "Set1"
#' 
#'   nMaxColors <- RColorBrewer::brewer.pal.info[base, 'maxcolors']
#'   if( is.null(n))
#'     n <- nMaxColors
#'   
#'   limit <- nMaxColors*(nMaxColors-1)/2
#'   if (n > limit){
#'     stop('Number of colors exceed limit of ', limit, ' colors per palette.')
#'   }
#'   
#'   if(n > nMaxColors){
#'     palette <- RColorBrewer::brewer.pal(nMaxColors, base)
#'     allComb <- combn(palette, 2)
#'     
#'     for (i in 1:(n-nMaxColors))
#'       palette <- c(palette, grDevices::colorRampPalette(allComb[,i])(3)[2])
#'     
#'   } else {
#'     palette <- RColorBrewer::brewer.pal(nMaxColors, base)[1:n]
#'   }
#'   palette
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' #' @title Builds a complete color palette for the conditions given in argument
#' #' 
#' #' @description xxxx
#' #' 
#' #' @param conds The extended vector of samples conditions
#' #' 
#' #' @param palette A vector of HEX color code that form the basis palette from which to build
#' #' the complete color vector for the conditions.
#' #' 
#' #' @return A vector composed of HEX color code for the conditions
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' #' conditions <- colData(Exp1_R25_pept)$Condition
#' #' GetColorsForConditions(conditions, ExtendPalette(2))
#' #' 
#' #' @export
#' #' 
#' #' @importFrom RColorBrewer brewer.pal
#' #' 
#' #' @return NA
#' #' 
#' GetColorsForConditions <- function(conds, palette=NULL){
#' 
#'   if(missing(conds))
#'     stop("'conds' is required")
#'   
#' 
#'   if (!is.null(palette) && length(unique(conds)) != length(palette))
#'     stop('The length of `conds` must be equal to the length of `base_palette`.')
#'   
#'   
#'   if (is.null(palette))
#'     palette <- ExtendPalette(length(unique(conds)))
#' 
#'   
#'   myColors <- NULL
#'   for (i in 1:length(conds))
#'     myColors[i] <- palette[which(conds[i] == unique(conds))]
#' 
#'   return(myColors)
#'   
#' }  
#' 
#' 
#' 
#' #' @title
#' #' xxx
#' #' 
#' #' @description
#' #' Builds a vector of colors (one per condition) to be used as a basis palette
#' #' 
#' #' @param palette.name The name of the palette from the package `RColorBrewer`
#' #' 
#' #' @param conditions xxxx
#' #' 
#' #' @importFrom RColorBrewer brewer.pal
#' #' 
#' #' @export
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' #' conditions <- colData(Exp1_R25_pept)$Condition
#' #' Base_Palette(palette.name = 'Dark2', conditions)
#' #'
#' #' @return NA
#' #' 
#' Base_Palette <- function(palette.name = 'Dark2', conditions){
#'   
#'   nbConds <- length(unique(conditions))
#'   nbColors <- max(3, nbConds)
#'   basePalette <- RColorBrewer::brewer.pal(nbColors, palette.name)[1:nbConds]
#'   return(basePalette)
#' }