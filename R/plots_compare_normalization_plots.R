#' Plot to compare the quantitative proteomics data before and after 
#' normalization using the library \code{highcharter}
#' 
#' @title Builds a plot from a numeric matrix.
#' @param qDataBefore A numeric matrix that contains quantitative data before 
#' normalization.
#' @param qDataAfter A numeric matrix that contains quantitative data after 
#' normalization.
#' @param condsForLegend A vector of the conditions (one condition 
#' per sample).
#' @param indData2Show A vector of the indices of the columns to show in 
#' the plot. The indices are those of indices of 
#' the columns int the data.frame qDataBefore.
#' @param palette xxx
#' @return A plot
#' @author Samuel Wieczorek, Enora Fremy
#' @examples
#' library(highcharter)
#' library(DAPAR)
#' # source('~/TELETRAVAIL/github_DAPARforFeatures/DAPAR2/R/normalize.R')
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' conds <- colData(obj)[["Condition"]]
#' objBefore <- obj[['original_log']]
#' qDataBefore <- assay(objBefore)
#' objAfter <- wrapper.normalizeD(objBefore,"QuantileCentering", conds, "within conditions")
#' qDataAfter <- assay(objAfter)
#' compareNormalizationD_HC(qDataBefore, qDataAfter, condsForLegend=conds)
#' @importFrom RColorBrewer brewer.pal
#' @import highcharter
#' @export
compareNormalizationD_HC <- function(qDataBefore,
                                     qDataAfter,
                                     condsForLegend=NULL,
                                     indData2Show=NULL,
                                     palette = NULL){
  
  if (is.null(condsForLegend)) return(NULL)
  if (is.null(indData2Show)) {indData2Show <- c(1:ncol(qDataAfter)) }
  
  
  if (is.null(palette)){
    tmp <- RColorBrewer::brewer.pal(length(unique(condsForLegend)),"Dark2")[1:length(unique(condsForLegend))]
    for (i in 1:ncol(qDataBefore)){
      palette[i] <- tmp[ which(condsForLegend[i] == unique(condsForLegend))]
    }
  }else{
    if (length(palette) != ncol(qDataBefore)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
  }
  
  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  ##Colors definition
  legendColor <- unique(palette)
  txtLegend <- unique(condsForLegend)
  
  
  series <- list()
  for (i in 1:length(indData2Show)){
    tmp <- list(name=condsForLegend[i], data =list_parse(data.frame(x=x[,indData2Show[i]],
                                                                    y=y[,indData2Show[i]])))
    series[[i]] <- tmp
  }
  
  h1 <-  highchart() %>% 
    dapar_hc_chart( chartType = "scatter") %>%
    hc_add_series_list(series) %>%
    hc_tooltip(enabled= "false" ) %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  
  return(h1)
  
}

