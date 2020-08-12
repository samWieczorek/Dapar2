#' The filtering functions of DAPAR have not been all moved to DAPAR2 as we now use the QFeatures package
#' which provides some filtering functions, especially on features that are present in the rowData
#' of the datasets.
#'  The filtering functions on numerical values are deleted because the same functions exist in QFeatures.
#'  For the missing values filtering on conditions, we do not use DAPAR functions anymore. Instead, we use
#'  the numerical filtering functions in QFeatures. To do so, it is necessary to build some rowdata for the
#'  SummarizedExperiment (not necessary stored in the object) which count the number of missing values w.r.t.
#'  the type of filtering: whole matrix, at least one value per condition, etc...
#'  
#'  
#'  






#' getPourcentageOfMV <- function(obj)
#' Voir nNA(object, i) de la classe QFeatures






#' @title Barplot of proportion of contaminants and reverse
#' 
#' @param lDataset The total length (number of rows) of the dataset
#' 
#' @param nBoth The number of both contaminants and reverse identified in the dataset.
#' 
#' @param nCont The number of contaminants identified in the dataset.
#' 
#' @param nRev The number of reverse entities identified in the dataset.
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





#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @description Returns the \code{SummarizedExperiment} object with a extra column in \code{rowData()}.
#' The extra column indicates feature(s) to delete in 1 and 0 otherwise.
#' The user chooses the threshold, either the percentage of NAS per row (default)
#' or the number of samples containing NAs per row allowed. Then, 
#' the filter tags the lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#'
#' The different methods are :
#' - "WholeMatrix": given a threshold \code{th}, only the lines that contain
#'   at least \code{th} values are kept.
#' - "AllCond": given a threshold \code{th}, only the lines which contain
#'   at least \code{th} values for each of the conditions are kept.
#' - "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#'   at least \code{th} values, and for at least one condition, are kept.
#'
#' @title An object of class \code{SummarizedExperiment} where the rowData gets an extra column 
#' containing the information of the fitration.
#'
#' @param object An object of class \code{QFeatures}
#' 
#' @param sampleTab \code{colData()} of obj
#' 
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' 
#' @param th Either a numeric between 0 and 1 where only the lines which contain
#' at least \code{th}% of non NA values are kept.
#' Or a integer between 0 and maximum number of samples for 'WholeMatrix' and 
#' between 0 and maximum number of replicate for "AllCond" and "AtLeastOneCond", where
#' only the lines which contain at least \code{th} values are kept.
#' 
#' @param percent TRUE by default. When FALSE, use the number of samples
#' 
#' @param newColName Name of the new column containing the filtration information
#' 
#' @return The object of class \code{SummarizedExperiment} with extra column in rowData
#' indicating 1 for the lines to remove, else 0.
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[100:120,]
#' res <- MVrowsTagToOne(obj, 2, type = "AtLeastOneCond", th=3, percent=FALSE, newColName="tagNA")
#' names(rowData(res))
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' 
MVrowsTagToOne <- function(object, type, th=0, percent=TRUE, newColName="tagNA") {
  
  if (is.null(object)) { return(NULL) }
  
  i <- length(experiments(object))
  
  ## Create a fictive coluln for each assays otherwise the filter method
  ## of QFeatures will truncate the object
  for (k in 1:length(experiments(object)))
    rowData(object[[k]])[newColName] <- rep(0, nrow(object[[k]]))
  
  sampleTab <- colData(object)
  
  # # check Type param
  # paramtype<-c("None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond") 
  # if (sum(is.na(match(type, paramtype)==TRUE))>0)
  #   stop("Param type is not correct.")
  
  # # check Threshold/Percent param
  # if (!is.numeric(th)) stop("th must be numeric.")
  # if (th < 0) {
  #   warning("Param th can't be inferior to zero.")
  #   th<-0
  #   cat("th set to 0.\n")
  # }
  # if (th==0) { cat("All row with NAs are kept.\n") }
  
  # if (isTRUE(percent)) {
  #   if (th>1) {
  #     th<-1
  #     warning("When percent=T, th can't be superior to one.")
  #     cat("th set to 1.\n")
  #   }
  #   if (th != 0) { cat("Row(s) containing",round(th, digits = 2)*100,"% NAs or more are removed.\n") 
  #   }
  # }
  # else { # !isTRUE(percent)
  #   if (th!=0) {
  #     if (th%%1 != 0) {
  #       stop("Param th have to be an integer.")
  #     }
  #     if (th > nrow(sampleTab)) {
  #       warning("When percent=F, param th can't be superior to the number of samples.")
  #       th<-nrow(sampleTab)
  #       cat("th set to the number of samples.\n")
  #     }
  #     cat("Rows are removed when at least",th,"sample(s) contain NAs.\n")
  #   }
  # }
  
  # Filtration
  keepThat <- NULL
  if (is.null(metadata(object[[i]])$OriginOfValues)) {
    data <- as.data.frame(assay(object[[i]]))
  } else {
    data <- dplyr::select(rowData(object[[i]]),metadata(object[[i]])$OriginOfValues)
  }
  
  if (type == "None") {
    keepThat <- seq(1:nrow(data))
  } else if (type == "EmptyLines") {
    keepThat <- which(apply(!is.MV(data), 1, sum) >= 1)
  } else if (type == "WholeMatrix") {
    if (isTRUE(percent)) {
      keepThat <- which(rowSums(!is.MV(data))/ncol(data) >= th) 
    } else {
      keepThat <- which(apply(!is.MV(data), 1, sum) >= th)
    }
  } else if (type == "AtLeastOneCond" || type == "AllCond") {
    
    if (is.null(sampleTab)) { return(NULL) }
    
    conditions <- unique(sampleTab[['Condition']])
    nbCond <- length(conditions)
    s <- matrix(rep(0, nrow(data)*nbCond),nrow=nrow(data), ncol=nbCond)
    
    if (isTRUE(percent)) {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        s[,c] <- (rowSums(!is.MV(data[,ind]))/length(ind)) >= th
      }
    } else {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        if (length(ind) == 1){
          s[,c] <- !is.MV(data[,ind]) >= th 
        }
        else {
          s[,c] <- (apply(!is.MV(data[,ind]), 1, sum)) >= th
        }
      }
    }
    
    if (type == "AllCond") {
      keepThat <- which(rowSums(s) == nbCond)
    }
    else if (type == "AtLeastOneCond") {
      keepThat <- which(rowSums(s) >= 1)
    }
  }
  
    rowData(object[[i]])[[newColName]][-keepThat] <- 1 
  
  return(object)
}



#' Returns the \code{SummarizedExperiment} object without the \code{rowData()} extra column
#' from MVrowsTagToOne
#'
#' @title Restore the rowData() before MVrowsTagToOne
#' 
#' @param object An object of class \code{QFeatures}
#' 
#' @param colToRemove Column to remove
#' 
#' @return The \code{SummarizedExperiment} object without the rowData colToRemove column
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[[2]]
#' sampleTab <- colData(Exp1_R25_pept)
#' obj <- MVrowsTagToOne(obj=object, i=2, sampleTab, newColName="LinesKept", type="AllCond", th=2)
#' obj <- removeAdditionalCol(obj, i=2, colToRemove="LinesKept" )
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' 
removeAdditionalCol <- function(object, colToRemove=NULL) {
  
  if(is.null(object)) { return(NULL) }
  if(is.null(colToRemove)){return(NULL)} 
  
  
  for (k in 1:length(experiments(object)))
    {
    if( is.na(match(colToRemove,names(rowData(object[[k]])))) ) {
      print(paste0("Warning: ",colToRemove," isn't a column name of rowData(obj)"))
    }
    else 
      rowData(object[[k]]) <- rowData(object[[k]])[,-which(colnames(rowData(object[[k]]))==colToRemove)]
  }
  
  return(object)
  
}