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



#' Returns the \code{SummarizedExperiment} object with a extra column in \code{rowData()}.
#' The extra column indicates feature(s) to delete in 1 and 0 otherwise.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter tags the lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' @param obj An object of class \code{SummarizedExperiment} containing
#' quantitative data.
#' @param newColName Name of the new column in rowData()
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' @param th An integer value of the threshold
#' @return The object of class \code{SummarizedExperiment} with extra column in rowData
#' indicating 0 for the lines to keep, else 1.
#' @author Florence Combes, Samuel Wieczorek, Enora Fremy
#' @examples
#' library(DAPAR2)
#' library(Features)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[[2]]
#' sampleTab <- colData(Exp1_R25_pept)
#' obj<-MVindicesToZero(obj=object, sampleTab, newColName="LinesKept", type="AllCond", th=2)
#' head(rowData(obj))
#' @export
MVindicesToZero <- function(obj, sampleTab, newColName, type, th = 0) {
  
  #Check parameters
  paramtype<-c("None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond") 
  if (sum(is.na(match(type, paramtype)==TRUE))>0)
    stop("Param type is not correct.")
  
  
  paramth<-c(seq(0, nrow(sampleTab), 1))
  if (sum(is.na(match(th, paramth)==TRUE))>0){
    warning("Param th is not correct.")
    return (NULL)
  }
  
  keepThat <- NULL
  if (is.null(metadata(obj)$OriginOfValues)){
    data <- as.data.frame(assay(obj))
  } else {
    data <- dplyr::select(rowData(obj),metadata(obj)$OriginOfValues)
  }
  
  if (type == "None"){
    keepThat <- seq(1:nrow(data))
  } else if (type == "EmptyLines"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= 1)
  } else if (type == "WholeMatrix"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= th)
  } else if (type == "AtLeastOneCond" || type == "AllCond"){
    
    conditions <- unique(sampleTab[['Condition']])
    nbCond <- length(conditions)
    keepThat <- NULL
    s <- matrix(rep(0, nrow(data)*nbCond),nrow=nrow(data), 
                ncol=nbCond)
    
    for (c in 1:nbCond){
      ind <- which(sampleTab[['Condition']] == conditions[c])
      if (length(ind) == 1){
        s[,c] <- (!is.MV(data[,ind]) >= th)}
      else {
        s[,c] <- (apply(!is.MV(data[,ind]), 1, sum) >= th)
      }
    }
    
    if (type == "AllCond") {
      keepThat <- which(rowSums(s) == nbCond)
    }
    else if (type == "AtLeastOneCond") {
      keepThat <- which(rowSums(s) >= 1)
    }
  }
  
  rowData(obj)[[newColName]]<-1
  rowData(obj)[[newColName]][keepThat]<-0
  
  return(obj)
}


#' Returns the \code{SummarizedExperiment} object without the \code{rowData()} extra column
#' from MVindicesToZero
#'  
#' @title Restore the rowData() before MVindicesToZero
#' @param obj An object of class \code{SummarizedExperiment}
#' @param newColName Column to remove
#' @return The \code{SummarizedExperiment} object without the rowData newColName column
#' @author Enora Fremy
#' @examples
#' library(DAPAR2)
#' library(Features)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[[2]]
#' sampleTab <- colData(Exp1_R25_pept)
#' obj<-MVindicesToZero(obj=object, sampleTab, newColName="LinesKept", type="AllCond", th=2)
#' removeAdditionalCol(obj,newColName="LinesKept" )
#' @export
removeAdditionalCol <- function(obj, newColName) {
  
  if( is.na(match(newColName,names(rowData(obj)))) ){ 
    print(paste0("Warning: ",newColName," isn't a column name of rowData(obj)"))
  }
  
  rowData(obj)[[newColName]] <- NULL
  
  return(obj)
}