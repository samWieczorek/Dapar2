#' Builds a densityplot of the CV of entities in numeric matrix.
#' The CV is calculated for each condition present in the dataset
#' (see the slot \code{'Condition'} in the \code{colData()} DataFrame)
#' 
#' @title Distribution of CV of entities
#' 
#' @description 
#' Builds a densityplot of the CV of entities in numeric matrix.
#' The CV is calculated for each condition present in the dataset
#' (see the slot \code{'Condition'} in the \code{colData()} DataFrame)
#' 
#' @param qData A numeric matrix that contains quantitative data.
#' 
#' @param conds A vector of the conditions (one condition per sample).
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
#' qData <- assay(Exp1_R25_pept[[2]])[1:10,]
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' CVDistD_HC(qData, conds)
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @import highcharter
#' 
#' @importFrom DT JS
#' 
#' @importFrom stats var
#' 
#' @export
#' 
CVDistD_HC <- function(qData, 
                       conds, 
                       palette = NULL){
  
  if (is.null(conds)) {return(NULL)}
  conditions <- unique(conds)
  n <- length(conditions)
  
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
    dapar_hc_chart(chartType = "spline", zoomType="x") %>%
    highcharter::hc_colors(unique(myColors)) %>%
    highcharter::hc_legend(enabled = TRUE) %>%
    highcharter::hc_xAxis(title = list(text = "CV(log(Intensity))")) %>%
    highcharter::hc_yAxis(title = list(text = "Density")) %>%
    highcharter::hc_tooltip(headerFormat= '',
                            pointFormat = "<b>{series.name}</b>: {point.y} ",
                            valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "logIntensity") %>%
    highcharter::hc_plotOptions(
      series=list(
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
      )
    )
  
  minX <- maxX <- 0
  maxY <- 0
  for (i in 1:n){
    if (length(which(conds == conditions[i])) > 1){
      t <- apply(qData[,which(conds == conditions[i])], 1, 
                 function(x) 100*stats::var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      tmp <- data.frame(x = density(t, na.rm = TRUE)$x,
                        y = density(t, na.rm = TRUE)$y)
      
      ymaxY <- max(maxY,tmp$y)
      xmaxY <- tmp$x[which(tmp$y==max(tmp$y))]
      minX <- min(minX, tmp$x)
      maxX <- max(maxX, 10*(xmaxY-minX))
      
      
      h1 <- h1 %>% hc_add_series(data=tmp, name=conditions[i]) }
  }
  
  h1 <- h1 %>%
    hc_chart(
      events = list(
        load = DT::JS(paste0("function(){
                         var chart = this;
                         this.xAxis[0].setExtremes(",minX,",",maxX, ");
                         this.showResetZoom();}"))
      )
    )
  
  return(h1)
  
}