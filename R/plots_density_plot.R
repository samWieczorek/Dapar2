#' @title Builds a densityplot of the quantitative data
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
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' densityPlotD_HC(qData, conds)
#' 
#' legend <- colData(Exp1_R25_pept)[["Sample.name"]]
#' pal <- ExtendPalette(2, 'Dark2')
#' densityPlotD_HC(qData, conds, legend, palette=pal)
#' 
#' @import highcharter
#' @importFrom stats density
#' 
#' @export
#' 
densityPlotD_HC <- function(qData, 
                            conds, 
                            legend=NULL, 
                            palette = NULL){
  
  if(missing(qData))
    stop("'qData' is missing.")
  
  if(missing(conds))
   stop("'conds' is missing.")
  
  if (is.null(legend))
    legend <- conds

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
    ) %>% hc_colors(myColors)
  
 
  
  for (i in 1:ncol(qData)){
    
    tmp <- data.frame(x = stats::density(qData[,i], na.rm = TRUE)$x, 
                      y = stats::density(qData[,i], na.rm = TRUE)$y)
    
    h1 <- h1 %>% hc_add_series(data=list_parse(tmp), name=legend[i]) 
  }

  return(h1)
}