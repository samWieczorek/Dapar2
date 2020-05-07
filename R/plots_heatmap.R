

#' Heatmap of the quantitative proteomic data of a \code{data.frame} object
#'
#' @title This function is a wrapper to \code{\link{heatmap.2}} that displays
#' a numeric matrix
#' 
#' @param qData A dataframe of numeric values
#' 
#' @param distance The distance used by the clustering algorithm to compute
#' the dendrogram. See \code{help(heatmap.2)}
#' 
#' @param cluster the clustering algorithm used to build the dendrogram.
#' See \code{help(heatmap.2)}
#' 
#' @param dendro A boolean to indicate if the dendrogram has to be displayed
#' 
#' @return A heatmap
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' ft <- Exp1_R25_pept
#' qData <- assay(ft[['original_log']])[1:1000,]
#' heatmapD(qData)
#' }
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom gplots heatmap.2
#' @importFrom stats dist hclust
#' @import graphics
#' 
#' @export
#' 
heatmapD <- function(qData, distance="euclidean", cluster="complete", dendro = FALSE){
  ##Check parameters
  # paramdist <- c("euclidean", "manhattan")
  # if (!(distance %in% paramdist)){
  #     stop("Param distance is not correct.")
  #     return (NULL)
  # }
  #
  # paramcluster <- c("ward.D", "average")
  # if (!(cluster %in%  paramcluster)){
  #     stop("Param clustering is not correct.")
  #     return (NULL)
  # }
  
  
  # if (isTRUE(dendro) && getNumberOfEmptyLines(qData) != 0)  {
  #     stop("Your dataset contains empty lines: the dendrogram cannot
  # be computed.
  #         Please filter or impute missing values before.")
  #     return (NULL)
  # }
  # else {
  .data <- matrix(qData,
                  ncol = ncol(qData),
                  byrow = FALSE,
                  dimnames = list(rownames(qData), colnames(qData))
  )
  colors = c(seq(-3, -2, length=100),
             seq(-2, 0.5, length=100),
             seq(0.5, 6, length=100))
  heatmap.color <- grDevices::colorRampPalette(c("green", "red"))(n = 1000)
  
  
  if (dendro){ .dendro = "row"} else {.dendro = "none"}
  p <- gplots::heatmap.2(
    x=t(.data),
    distfun = function(x) {
      x[is.na(x)] <- -1e5
      dist(x, method=distance)
    },
    hclustfun = function(x) {
      x[is.na(x)] <- -1e5
      hclust(x, method=cluster)
    },
    dendrogram =.dendro,
    Rowv=TRUE,
    col=heatmap.color ,
    density.info='none',
    key=TRUE,
    trace="none",
    scale="none",
    #srtCol=45,
    labCol="",
    margins=c(4,12),
    cexRow=1.5,
    keysize = 1.5,
    lhei = c(1.5, 9),
    lwid = c(1.5, 4),
    lmat = rbind(4:3, 2:1)
    
  )
  #    }
}

