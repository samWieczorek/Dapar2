#' Densityplot of quantitative proteomics data over samples. 
#' 
#' @title Builds a densityplot from a dataframe
#' 
#' @param qData numeric matrix
#' 
#' @param conds DataFrame
#' 
#' @param legend A vector of the conditions (one condition per sample).
#' 
#' @param palette A vector of HEX Code colors (one color per condition).
#' 
#' @return A density plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' legend <- colData(Exp1_R25_pept)[["Sample.name"]]
#' palette <- c("#fc03ec","#94fc03")
#' densityPlotD_HC(qData, conds, legend, palette)
#' 
#' @import highcharter
#' @importFrom stats density
#' 
#' @export
densityPlotD_HC <- function(qData, conds, legend=NULL, palette = NULL){
  
 if(is.null(conds))
   stop("'conds' is missing.")
  
  palette <- BuildPalette(conds, palette)
  
  h1 <-  highcharter::highchart() %>% 
    hc_title(text = "Density plot") %>% 
    dapar_hc_chart(chartType = "spline", zoomType="x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_xAxis(title = list(text = "log(Intensity)")) %>%
    hc_yAxis(title = list(text = "Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "densityplot") %>%
    hc_plotOptions(
      series=list(
        animation=list(
          duration = 100
        ),
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
      )
    )
  
  if (!is.null(palette)) {
    if (length(palette) != ncol(qData)){
      warning("The color palette has not the same dimension as the number of samples")
      return(NULL)
    }
    h1 <- h1 %>% hc_colors(palette)
  }
  
  if (is.null(legend)) {
    legend <- paste0("series", 1:ncol(qData))
  }
  
  for (i in 1:ncol(qData)){
    
    tmp <- data.frame(x = stats::density(qData[,i], na.rm = TRUE)$x, 
                      y = stats::density(qData[,i], na.rm = TRUE)$y)
    
    h1 <- h1 %>% hc_add_series(data=list_parse(tmp), name=legend[i]) 
    
  }
  
  
  return(h1)
  
}