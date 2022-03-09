#' @param qData numeric matrix
#' 
#' @param conds A `character()` of condition name for each sample. The 
#' length of 'conds' must be equal to the number of columns of 'qData'.
#' 
#' @param legend A vector of the conditions (one condition per sample).
#' 
#' @param pal.name A `character(1)` which is the name of a palette in
#' the package [RColorBrewer].
#' 
#' @import highcharter
#' @importFrom stats density
#' 
#' @export
#' 
#' @rdname descriptive-statistics
#' 
densityPlot <- function(object, 
                        conds, 
                        pal.name = NULL){
  
  if(missing(qData))
    stop("'qData' is missing.")
  
  if(missing(conds))
   stop("'conds' is missing.")
  
  if (length(conds) != ncol(qData))
    stop("qData and conds must have the same number of samples.")
  
  if (is.null(pal.name)){
    warning("Color palette set to default.")
    myColors <- SampleColors(conds)
  } else
    myColors <- SampleColors(conds, pal.name)

 
  h1 <-  highcharter::highchart() %>% 
    hc_title(text = "Density plot") %>% 
    customChart(chartType = "spline", zoomType="x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_xAxis(title = list(text = "log(Intensity)")) %>%
    hc_yAxis(title = list(text = "Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    customExportMenu(fname = "densityplot") %>%
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
  
 
  
  for (i in seq_len(ncol(qData))){
    
    tmp <- data.frame(x = stats::density(qData[,i], na.rm = TRUE)$x, 
                      y = stats::density(qData[,i], na.rm = TRUE)$y)
    
    h1 <- h1 %>% 
      hc_add_series(data=list_parse(tmp), 
                    name = colnames(qData)[i]) 
  }

  h1
}