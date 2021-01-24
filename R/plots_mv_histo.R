
#' This method plots a histogram of missing values.
#' 
#' @title Histogram of missing values
#' 
#' @param qData A dataframe that contains quantitative data.
#' 
#' @param conds A vector of the conditions (one condition per sample).
#' 
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' 
#' @param palette xxx
#' 
#' @return A histogram
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' conds <- colData(Exp1_R25_pept)
#' mvHisto_HC(qData, conds, showValues=TRUE)
#' 
#' pal <- ExtendPalette(2, 'Dark2')
#' mvHisto_HC(qData, conds, showValues=TRUE, palette = pal)
#' 
#' 
#' @export
#' 
#' @import highcharter
#' 
mvHisto_HC <- function(qData, 
                       conds, 
                       showValues = FALSE, 
                       palette = NULL){
  
  myColors <- NULL
  if (is.null(palette)){
    warning("Color palette set to default.")
    myColors <-  GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
  } else {
    if (length(palette) != length(unique(conds))){
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- GetColorsForConditions(conds, ExtendPalette(length(unique(conds))))
    } else 
      myColors <- palette
    myColors <- GetColorsForConditions(conds, palette)
  }
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  
  df <- data.frame(NbNAPerCol)
  names(df) <- 'y'
  
  
  h1 <-  highchart() %>%
    dapar_hc_chart(chartType = "column") %>%
    hc_title(text = "#NA by replicate") %>%
    hc_add_series(df,type="column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds[['Condition']], title = list(text = "Replicates")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y}")
  
  
  return(h1)
 
}

