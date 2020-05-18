#' The filtering functions of DAPAR have not been all moved to DAPAR2 as we now use the Features package
#' which provides some filtering functions, especially on features that are present in the rowData
#' of the datasets.
#'  The filtering functions on numerical values are deleted because the same functions exist in Features
#'  For the missing values filtering on conditions, we do not use DAPAR functions anymore. Instead, we use
#'  the numerical filtering functions in Features. To do so, it is necessary to build some rowdata for the
#'  SummarizedExperiment (not necessary stored in the object) which count the number of missing values w.r.t.
#'  the type of filtering: whole matrix, at least one value per condition, etc...
#'  
#'  
#'  






#' getPourcentageOfMV <- function(obj)
#' Voir nNA(object, i) de la classe Features






#' @title Barplot of proportion of contaminants and reverse
#' 
#' @param nBoth The number of both contaminants and reverse identified in the dataset.
#' 
#' @param nCont The number of contaminants identified in the dataset.
#' 
#' @param nRev The number of reverse entities identified in the dataset.
#' 
#' @param lDataset The total length (number of rows) of the dataset
#' 
#' @return A barplot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' proportionConRev_HC(100, 10, 20)
#' 
#' @export
#' 
#' @import highcharter
#' 
proportionConRev_HC <- function(lDataset, nBoth = 0, nCont=0, nRev=0){
  
  if(missing(lDataset))
    stop("'lDataset' is missing.")
  
  
  total <- nBoth + nCont + nRev + lDataset
  pctGood <- 100 * round(lDataset/total,  digits=4)
  pctBoth <- 100 * round(nBoth/total,  digits=4)
  pctContaminants <- 100 * round(nCont/total,  digits=4)
  pctReverse <- 100 * round(nRev/total,  digits=4)
  
  counts <- c(lDataset, nCont, nRev, nBoth)
  slices <- c(pctGood, pctContaminants, pctReverse ,pctBoth) 
  lbls <- c("Quantitative data", "Contaminants", "Reverse", "Both contaminants & Reverse")
  lbls <- paste(lbls, " (", counts, " lines)", sep="") 
  
  mydata <- data.frame(test=c(pctGood, pctContaminants, pctReverse ,pctBoth))
  
  highchart() %>% 
    dapar_hc_chart(chartType = "bar") %>% 
    hc_yAxis(title = list(text = "Pourcentage")) %>% 
    hc_xAxis(categories=lbls) %>% 
    hc_legend(enabled = FALSE) %>%
    hc_plotOptions(column = list(
      dataLabels = list(enabled = TRUE),
      stacking = "normal",
      enableMouseTracking = FALSE)
    ) %>% 
    hc_add_series(data  = mydata$test,
                  dataLabels = list(enabled = TRUE, format='{point.y}%'),
                  colorByPoint = TRUE) %>%
    dapar_hc_ExportMenu(filename = "contaminants")
  
  
}



#' Returns the indices of the lines of \code{exprs()} table to delete w.r.t. 
#' the conditions on the number of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "wholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "allCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "atLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' @param obj An object of class \code{SummarizedExperiment} containing
#' quantitative data.
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "wholeMatrix", "allCond", "atLeastOneCond"
#' @param th An integer value of the threshold
#' @return An vector of indices that correspond to the lines to keep.
#' @author Florence Combes, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mvFilterGetIndices(Exp1_R25_pept, "wholeMatrix", 2)
#' @export
mvFilterGetIndices <- function(obj, type, th = 0)
{
  #Check parameters
  paramtype<-c("None", "EmptyLines", "wholeMatrix", "allCond", "atLeastOneCond") 
  if (sum(is.na(match(type, paramtype)==TRUE))>0)
    stop("Param type is not correct.")

  
  paramth<-c(seq(0, nrow(Biobase::pData(obj)), 1))
  if (sum(is.na(match(th, paramth)==TRUE))>0){
    warning("Param th is not correct.")
    return (NULL)
  }
  
  keepThat <- NULL
  if (is.null(obj@experimentData@other$OriginOfValues)){
    data <- Biobase::exprs(obj)
  } else {
    data <- dplyr::select(Biobase::fData(obj),obj@experimentData@other$OriginOfValues)
  }
  
  if (type == "None"){
    keepThat <- seq(1:nrow(data))
  } else if (type == "EmptyLines"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= 1)
  } else if (type == "wholeMatrix"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= th)
  } else if (type == "atLeastOneCond" || type == "allCond"){
    
    conditions <- unique(Biobase::pData(obj)$Condition)
    nbCond <- length(conditions)
    keepThat <- NULL
    s <- matrix(rep(0, nrow(data)*nbCond),nrow=nrow(data), 
                ncol=nbCond)
    
    for (c in 1:nbCond){
      ind <- which(Biobase::pData(obj)$Condition == conditions[c])
      if (length(ind) == 1){
        s[,c] <- (!is.MV(data[,ind]) >= th)}
      else {
        s[,c] <- (apply(!is.MV(data[,ind]), 1, sum) >= th)
      }
    }
    
    
    if (type == "allCond") {
      keepThat <- which(rowSums(s) == nbCond)
    }
    else if (type == "atLeastOneCond") {
      keepThat <- which(rowSums(s) >= 1)
    }
  }
  return(keepThat)
}
