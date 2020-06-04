
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
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' mvHisto_HC(qData, conds, showValues=TRUE)
#' 
#' @export
#' 
#' @import highcharter
#' 
mvHisto_HC <- function(qData, conds, showValues = FALSE, palette = NULL){
  
  palette <- BuildPalette(conds, palette)
 
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  
  df <- data.frame(NbNAPerCol)
  names(df) <- 'y'
  
  
  h1 <-  highchart() %>%
    dapar_hc_chart(chartType = "column") %>%
    hc_title(text = "#NA by replicate") %>%
    hc_add_series(df,type="column", colorByPoint = TRUE) %>%
    hc_colors(palette) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds, title = list(text = "Replicates")) %>%
    dapar_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y}")
  
  
  return(h1)
 
}

