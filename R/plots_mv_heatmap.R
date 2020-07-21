#' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class \code{QFeatures} and 
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to vizualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' @title Heatmap of missing values
#' @param qData A matrix that contains quantitative data.
#' @param conds A vector of the conditions (one condition per sample).
#' @return A heatmap
#' @author Samuel Wieczorek, Thomas Burger, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' mvImage(qData, conds)
#' @export
#' @importFrom stats setNames
mvImage <- function(qData, conds){
  
  ### build indices of conditions
  indCond <- list()
  ConditionNames <- unique(conds)
  for (i in ConditionNames) {
    indCond <- append(indCond, list(which(i == conds)))
  }
  indCond <- stats::setNames(indCond, as.list(c("cond1", "cond2")))
  
  nNA1 = apply(as.matrix(qData[,indCond$cond1]), 1, function(x) sum(is.na(x)))
  nNA2 = apply(as.matrix(qData[,indCond$cond2]), 1, function(x) sum(is.na(x)))
  o <- order(((nNA1 +1)^2) / (nNA2 +1))
  exprso <- qData[o,]
  
  for (i in 1:nrow(exprso)){
    k <- order(exprso[i,indCond$cond1])
    exprso[i,rev(indCond$cond1)] <- exprso[i, k]
    .temp <- mean(exprso[i,rev(indCond$cond1)], na.rm = TRUE)
    exprso[i,which(!is.na(exprso[i,indCond$cond1]))] <- .temp
    
    k <- order(exprso[i,indCond$cond2])
    exprso[i,indCond$cond2] <- exprso[i, k+length(indCond$cond1)]
    .temp <- mean(exprso[i,indCond$cond2], na.rm = TRUE)
    exprso[i,length(indCond$cond1) + 
             which(!is.na(exprso[i,indCond$cond2]))] <- .temp
  }
  
  
  heatmap.DAPAR(exprso,
                col = colorRampPalette(c("yellow", "red"))(100),
                key=TRUE,
                srtCol= 0,
                labCol=conds,
                #ylab = "Peptides / proteins",
                main = "MEC heatmap"
  )
  
  #heatmap_HC(exprso,col = colfunc(100),labCol=conds)
  
  
}