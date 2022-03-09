#' @importFrom RColorBrewer brewer.pal
#' @import highcharter
#' @importFrom DT JS
#' @importFrom stats var
#' @importFrom SummarizedExperiment assay
#' 
#' @export
#' 
#' @rdname descriptive-statistics
#' 
CVDist <- function(object,
                   conds = NULL, 
                   pal.name = NULL){
  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  qData <- assay(object)
  
  
  if (is.null(conds))
    stop("'conds' is NULL")
  else if (length(conds) != ncol(qData))
    stop("'conds' must have the same length as the number of samples in 
         the dataset.")

  u_conds <- unique(conds)
  
  myColors <- SampleColors(u_conds)

  h1 <-  highcharter::highchart() %>% 
    customChart(chartType = "spline", zoomType="x") %>%
    highcharter::hc_colors(myColors) %>%
    highcharter::hc_legend(enabled = TRUE,
                           categories = u_conds) %>%
    highcharter::hc_xAxis(title = list(text = "CV(log(Intensity))")) %>%
    highcharter::hc_yAxis(title = list(text = "Density")) %>%
    highcharter::hc_tooltip(headerFormat= '',
                            pointFormat = "<b>{series.name}</b>: {point.y} ",
                            valueDecimals = 2) %>%
    customExportMenu(fname = "logIntensity") %>%
    highcharter::hc_plotOptions(
      series=list(
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
      )
    )
  
  minX <- maxX <- 0
  maxY <- 0
  for (i in seq_len(length(u_conds))){
    if (length(which(conds == u_conds[i])) > 1){
      t <- apply(qData[,which(conds == u_conds[i])], 1, 
                 function(x) 100*stats::var(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      tmp <- data.frame(x = density(t, na.rm = TRUE)$x,
                        y = density(t, na.rm = TRUE)$y)
      
      ymaxY <- max(maxY,tmp$y)
      xmaxY <- tmp$x[which(tmp$y == max(tmp$y))]
      minX <- min(minX, tmp$x)
      maxX <- max(maxX, 10*(xmaxY-minX))
      
      h1 <- h1 %>% hc_add_series(data=tmp, name=u_conds[i]) }
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