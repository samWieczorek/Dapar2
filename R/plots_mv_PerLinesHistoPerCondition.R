#' This method plots a bar plot which represents the distribution of the 
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#' 
#' @title Bar plot of missing values per lines and per condition
#' 
#' @param qData A dataframe that contains quantitative data.
#' 
#' @param samplesData A dataframe where lines correspond to samples and 
#' columns to the meta-data for those samples.
#' 
#' @param indLegend The indice of the column names of \code{colData()}
#' 
#' @param palette A vector of HEX Code colors (one color 
#' per condition).
#' 
#' @return A bar plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- SummarizedExperiment::assay(Exp1_R25_pept[[2]])
#' samplesData <- SummarizedExperiment::colData(Exp1_R25_pept)
#' mvPerLinesHistoPerCondition_HC(qData, samplesData)
#' 
#' @export
#' 
#' @import highcharter
#' 
mvPerLinesHistoPerCondition_HC <- function(qData, samplesData, indLegend="auto", palette=NULL){
  
  myColors <- NULL
  if (is.null(palette)){
    warning("Color palette set to default.")
    myColors <-   GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
  } else {
    if (length(palette) != length(unique(conds))){
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
    } else 
      myColors <- GetColorsForConditions(conds, palette)
  }
  
  if (identical(indLegend,"auto")) { indLegend <- c(2:length(colnames(samplesData)))}
  
  nbConditions <- length(unique(samplesData[["Condition"]]))
  
  ncolMatrix <- max(unlist(lapply(unique(samplesData[["Condition"]]), function(x){length(which(samplesData[["Condition"]]==x))})))
  m <- matrix(rep(0, nbConditions*(1+ncolMatrix)), 
              ncol = nbConditions, 
              dimnames=list(seq(0:(ncolMatrix)),unique(samplesData[["Condition"]])))
  
  for (i in unique(samplesData[["Condition"]])) {
    nSample <- length(which(samplesData[["Condition"]] == i))
    t <- NULL
    if (nSample == 1) {
      t <- table(as.integer(is.na(qData[,which(samplesData[["Condition"]] == i)])))
    } else {t <- table(rowSums(is.na(qData[,which(samplesData[["Condition"]] == i)])))}
    
    m[as.integer(names(t))+1,i] <- t
  }
  
  m <- as.data.frame(m)
  
  rownames(m) <- 0:(nrow(m)-1)
  
  h1 <-  highchart() %>% 
    hc_title(text = "#[lines] with X NA values (condition-wise)") %>% 
    dapar_hc_chart(chartType = "column") %>%
    hc_plotOptions( column = list(stacking = ""),
                    dataLabels = list(enabled = FALSE),
                    animation=list(duration = 100)) %>%
    hc_colors(unique(myColors)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(m), title = list(text = "#[NA values] per line (condition-wise)")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot_2") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y} ")
  
  for (i in 1:nbConditions){
    h1 <- h1 %>% hc_add_series(data=m[,unique(samplesData[["Condition"]])[i]]) }
  
  
  return(h1)
  
}