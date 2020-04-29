
#' Correlation matrix based on a data.frame object. Same as the 
#' function \link{corrMatrixD} but uses the package \code{highcharter}
#' 
#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#' @param object The result of the \code{cor} function.
#' @param names xxxxx
#' @param rate The rate parameter to control the exponential law for 
#' the gradient of colors
#' @return A colored correlation matrix
#' @author Samuel Wieczorek, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])[1:1000,]
#' names <- colData(Exp1_R25_pept)
#' res <- cor(qData,use = 'pairwise.complete.obs')
#' corrMatrixD_HC(res)
#' @importFrom dplyr tbl_df mutate left_join
#' @importFrom tidyr gather
#' @import highcharter
#' @importFrom DT JS
#' @importFrom tibble tibble
#' @export
corrMatrixD_HC <- function(qData, names = NULL, rate = 0.5) {
  
  df <- as.data.frame(qData)
  
  if (!is.null(names)){
      colnames(df) <- names
  }
  
  is.num <- sapply(df, is.numeric)
  df[is.num] <- lapply(df[is.num], round, 2)
  dist <- NULL
  
  x <- y <- names(df)
  
  df <- dplyr::tbl_df(cbind(x = y, df)) %>% 
    tidyr::gather(y, dist, -x) %>% 
    dplyr::mutate(x = as.character(x),
                  y = as.character(y)) %>% 
    dplyr::left_join(tibble::tibble(x = y,
                                    xid = seq(length(y)) - 1), by = "x") %>% 
    dplyr::left_join(tibble::tibble(y = y,
                                    yid = seq(length(y)) - 1), by = "y")
  
  ds <- df %>% 
    dplyr::select_("xid", "yid", "dist") %>% 
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