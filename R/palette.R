#' @title Palette for samples.
#'
#' @description
#' 
#' Builds a vector of `#conditions` colors for a set of samples. One color is
#' given for a given condition. This function extends a base palette from the package [RColorBrewer] to 'n' colors. 
#' The colors in the returned palette are always in the same order
#'
#' @param n The number of desired colors in the palette
#'
#' @param pal.name The name of the palette of the package [RColorBrewer] from which
#' the extended palette is built. Default value is 'Set1'. A `character(1)` which is the name of the palette from the package
#' [RColorBrewer] from which the colors are taken.
#' 
#' @param conds A `character()` of conditions. The length of the vector is
#' the number of samples in the dataset.
#'
#' @return A `HEX()` color codes.
#' 
#' @name palette
#' 
#' @examples
#' 
#' #-----------------------------------------------
#' # Builds a palette for a dataset with 3 conditions
#' # of 3 samples each.
#' #------------------------------------------------
#' 
#' conds <- c(rep('A', 3), rep('B', 3), rep('C', 3))
#' SampleColors(conds)
#' SampleColors(conds, 'Dark2')
#' 
#' 
#' #-----------------------------------------------
#' # Extend the default palette to 12 colors
#' #-----------------------------------------------
#' 
#' ExtendPalette(12)
NULL


#' @rdname palette
#' @export
#'
SampleColors <- function(conds, pal.name = 'Set1'){
if (missing(pal.name) || is.null(pal.name))
  pal.name <- "Set1"

  palette.ex <- NULL
  n.colors <- length(unique(conds))
  base.pal <- ExtendPalette(n.colors, pal.name)
  for (i in seq_len(n.colors)){
    ind <- which(conds == unique(conds)[i])
    palette.ex[ind] <- base.pal[i]
  }

  return(palette.ex)
}


#' @export
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
#' @importFrom utils combn
#' 
#' @rdname palette
#'
ExtendPalette <- function(n, pal.name = "Set1"){

  stopifnot(is.numeric(n))
  
  if( !(pal.name %in% rownames(brewer.pal.info)))
    pal.name <- "Set1"
  
  extended.pal <- NULL
  
  nMaxColors <- RColorBrewer::brewer.pal.info[pal.name, 'maxcolors']
  
  limit <- nMaxColors*(nMaxColors-1)/2
  if (n > limit){
    stop('Number of colors exceed limit of ', limit, ' colors per palette.')
  }

  if(n > nMaxColors){
    extended.pal <- RColorBrewer::brewer.pal(nMaxColors, pal.name)
    allComb <- combn(extended.pal, 2)

    for (i in 1:(n-nMaxColors))
      extended.pal <- c(extended.pal, 
                   grDevices::colorRampPalette(allComb[,i])(3)[2])

  } else {
    extended.pal <- RColorBrewer::brewer.pal(nMaxColors, pal.name)[1:n]
  }
  extended.pal
}
 