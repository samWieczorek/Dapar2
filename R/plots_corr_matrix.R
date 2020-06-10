

#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#' 
#' @param res xxx
#' 
#' @param names xxxxx
#' 
#' @param rate The rate parameter to control the exponential law for 
#' the gradient of colors
#' 
#' @return A colored correlation matrix
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[[2]]
#' corrMatrixD_HC(obj)
#' 
#' @importFrom dplyr mutate left_join
#' @importFrom tidyr gather
#' @import highcharter
#' @importFrom DT JS
#' @importFrom tibble tibble as_tibble
#' 
#' @export
#' 
corrMatrixD_HC <- function(obj, names = NULL, rate = 0.5) {
  
  res <- cor(SummarizedExperiment::assay(obj), use = 'pairwise.complete.obs')
  df <- as_tibble(res)
  
  if (!is.null(names)){
      colnames(df) <- names
  } 
  
  is.num <- sapply(df, is.numeric)
  df[is.num] <- lapply(df[is.num], round, 2)
  dist <- NULL
  
  x <- y <- names(df)
  
  df <- tibble::as_tibble(cbind(x = y, df)) %>% 
    tidyr::gather(y, dist, -x) %>% 
    dplyr::mutate(x = as.character(x),
                  y = as.character(y)) %>% 
    dplyr::left_join(tibble::tibble(x = y,
                                    xid = seq(length(y)) - 1), by = "x") %>% 
    dplyr::left_join(tibble::tibble(y = y,
                                    yid = seq(length(y)) - 1), by = "y")
  
  ds <- df %>% 
    dplyr::select("xid", "yid", "dist") %>% 
    highcharter::list_parse2()
  
  fntltp <- DT::JS("function(){
                  return this.series.xAxis.categories[this.point.x] + ' ~ ' +
                         this.series.yAxis.categories[this.point.y] + ': <b>' +
                         Highcharts.numberFormat(this.point.value, 2)+'</b>';
               ; }")
  cor_colr <- list( list(0, '#FF5733'),
                    list(0.5, '#F8F5F5'),
                    list(1, '#2E86C1')
  )
  highcharter::highchart() %>% 
    dapar_hc_chart(chartType = "heatmap") %>% 
    hc_xAxis(categories = y, title = NULL) %>% 
    hc_yAxis(categories = y, title = NULL) %>% 
    hc_add_series(data = ds) %>% 
    hc_plotOptions(
      series = list(
        boderWidth = 0,
        dataConditions = list(enabled = TRUE),
        dataLabels = list(enabled = TRUE)
      )) %>% 
    hc_tooltip(formatter = fntltp) %>% 
    hc_legend(align = "right", layout = "vertical",
              verticalAlign="middle") %>% 
    hc_colorAxis(  stops= cor_colr,min=rate,max=1) %>%
    dapar_hc_ExportMenu(filename = "corrMatrix")
}