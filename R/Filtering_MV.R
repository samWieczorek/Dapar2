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


#' @title Check the validity of the experimental design
#'
#' @description
#'
#' This manual page describes the computation of statistical test using [QFeatures] objects. In the following
#' functions, if `object` is of class `QFeatures`, and optional assay
#' index or name `i` can be specified to define the assay (by name of
#' index) on which to operate.
#'
#' The following functions are currently available:
#'
#' - `compute.t.test(xxxxx)` xxxxx.
#'
#' - `compute.group.t.test(xxxxx)` xxxxx.
#'
#' - `limma.complete.test(object, sampleTab)` uses the package Limma
#'
#'
#' @details xxx
#'
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' object <- addListAdjacencyMatrices(object, 3)
#' sTab <- colData(object)
#' gttest.se <- t_test_sam(object[[3]], sTab, FUN = 'limma.complete.test')
#'
#' object <- t_test_sam(object, 3, name = "ttestAssay", FUN = 'compute.t.test', contrast = 'OnevsOne')
#'
"filterFeatures_sam"

#' @param  object An object of class `QFeatures`.
#'
#' @param name A `character(1)` naming the new assay name. Defaults
#'     are `ttestAssay`.
#'
#' @param FUN xxx
#'
#' @param ... Additional parameters passed to inner functions.
#'
#' @export
#'
#' @rdname filterFeatures_sam
#'
#' @example 
#' library(QFeatures)
#' utils::data(Exp1_R25_prot, package='DAPARdata2')
#' object <- Exp1_R25_prot
#' filter <- VariableFilter(field='Reverse', value='+', condition='!=')
#' object <- filterFeatures_sam(object, 2, filter)
#' 
setMethod("filterFeatures_sam", "QFeatures",
          function(object, i, name = "filterAssay", filter) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]

 

            argg <- c(as.list(environment()))
            tmp <-  filterFeatures_sam(object[[i]], filter)

            
            object <- addAssay(object,
                               tmp,
                               name)
            addAssayLinkOneToOne(object, from = i, to = name)
          })



#' @importFrom AnnotationFilter condition, field
setMethod("filterFeatures_sam", "SummarizedExperiment",
          function(object, filter) {
            x <- rowData(object)
            if (field(filter) %in% names(x)){
              sel <- do.call(condition(filter),
                      list(x[, field(filter)], value(filter))
                      )
            } else {
              sel <- rep(FALSE, nrow(x))
            }
            object[sel, ]
          }
          )



#' @importFrom BiocGenerics do.call
filterFeaturesWithAnnotationFilter <- function(object, filter, na.rm, ...) {
  sel <- lapply(experiments(object),
                function(exp) {
                  x <- rowData(exp)
                  if (field(filter) %in% names(x))
                    do.call(condition(filter),
                            list(x[, field(filter)],
                                 value(filter)))
                  else
                    rep(FALSE, nrow(x))
                })
  sel <- lapply(sel, function(x) {
    x[is.na(x)] <- !na.rm
    x
  })
  if (not(filter)) sel <- lapply(sel, "!")
  object[sel, ]
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
#' @return The object of class \code{SummarizedExperiment} with extra column in rowData
#' indicating 1 for the lines to remove, else 0.
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' res <- MVrowsTagToOne(Exp1_R25_pept, type = "WholeMatrix", th=1, percent=FALSE)
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' 
MVrowsTagToOne <- function(object=NULL, type=NULL, th=0, percent=TRUE) {
  
  if (is.null(object) || is.null(type)) { return(NULL) }
  
  if (missing(object))
    stop("'object' is required")
  if (missing(type))
    stop("'type' is required")
  
  if (!(type %in% c('None', 'EmptyLines', 'WholeMatrix', 'AtLeastOneCond', 'AllCond')))
      stop("Available values for 'type' are not one of: 'None', 'EmptyLines', 'WholeMatrix', 'AtLeastOneCond', 'AllCond'")
      
  newColName <- "tagNA"
  i <- length(experiments(object))
  
  ## Create a fictive column for each assays otherwise the filter method
  ## of QFeatures will truncate the object
  for (k in 1:length(experiments(object)))
    rowData(object[[k]]) <- setNames(cbind(rowData(object[[k]]),tmp=rep(0, nrow(object[[k]]))),
                                      c(names(rowData(object[[k]])),newColName)
    )
   
  
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
  keepThat <- df <- NULL
  #if (is.null(metadata(object[[i]])$OriginOfValues)) {
  #  print('ok')
    df <- as.data.frame(assay(object[[i]]))
  #} else {
  #  print('pas ok')
  #  df <- dplyr::select(rowData(object[[i]]),metadata(object[[i]])$OriginOfValues)
  #}
  
  
  if (type == "None") {
    keepThat <- seq(1:nrow(df))
  } else if (type == "EmptyLines") {
    keepThat <- which(apply(!is.na(df), 1, sum) >= 1)
  } else if (type == "WholeMatrix") {
    if (isTRUE(percent)) {
      keepThat <- which(rowSums(!is.na(df))/ncol(df) >= th) 
    } else {
      keepThat <- which(apply(!is.na(df), 1, sum) >= th)
    }
  } else if (type == "AtLeastOneCond" || type == "AllCond") {
    
    if (is.null(sampleTab)) { return(NULL) }
    
    conditions <- unique(sampleTab[['Condition']])
    nbCond <- length(conditions)
    s <- matrix(rep(0, nrow(df)*nbCond),nrow=nrow(df), ncol=nbCond)
    
    if (isTRUE(percent)) {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        s[,c] <- (rowSums(!is.na(df[,ind]))/length(ind)) >= th
      }
    } else {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        if (length(ind) == 1){
          s[,c] <- !is.na(df[,ind]) >= th 
        }
        else {
          s[,c] <- (apply(!is.na(df[,ind]), 1, sum)) >= th
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
  
  
  rowData(object[[i]])[-keepThat,newColName] <- 1 
  
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
removeAdditionalCol <- function(object=NULL, colToRemove=NULL) {
  
  if(is.null(object)) { return(NULL) }
  if(is.null(colToRemove)){return(NULL)} 
  
  if (missing(object))
    stop("'object' is required")
  if (missing(colToRemove))
    stop("'colToRemove' is required")
  
  
  
  colToRemove <- "tagNA"
  
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



