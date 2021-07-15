

#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#' 
#' @param obj xxx
#' 
#' @param names A vector of strings which will be used as legend. The length
#' of 'names' must be equal to the number of samples in the dataset.
#' 
#' @param rate The rate parameter to control the exponential law for 
#' the gradient of colors
#' 
#' @param showValues xxx
#' 
#' @return A colored correlation matrix
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' corrMatrixD_HC(Exp1_R25_pept[[2]])
#' 
#' @importFrom dplyr mutate left_join
#' @importFrom tidyr gather
#' 
#' @import highcharter
#' 
#' @importFrom DT JS
#' 
#' @importFrom tibble tibble as_tibble
#' 
#' @importFrom stats cor
#' 
#' @export
#' 
corrMatrixD_HC <- function(obj.se, 
                           names = NULL, 
                           rate = 0.5, 
                           showValues=TRUE) {
  

    if (class(obj.se) != 'SummarizedExperiment'){
      warning("'obj.se' is not a SummarizedExperiment object.")
      return(NULL)
    }
  res <- cor(SummarizedExperiment::assay(obj.se),
             use = 'pairwise.complete.obs')

   #df <- as.data.frame(res)
   df <- tibble::as_tibble(res)

  if (!is.null(names)){
    if (length(names) != ncol(df)){
      warning("The length of 'names' must be equal to the number of samples
      in the dataset")
      return(NULL)
    }
    colnames(df) <- names
  }
      
  
  is.num <- sapply(df, is.numeric)
  df[is.num] <- lapply(df[is.num], round, 2)
  dist <- NULL
  
  x <- y <- names(df)

  df <- tibble::as_tibble(cbind(x = y, df)) %>% 
     #df <- dplyr::tbl_df(cbind(x = y, df)) %>% 
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
        dataLabels = list(enabled = showValues)
      )) %>% 
    hc_tooltip(formatter = fntltp) %>% 
    hc_legend(align = "right", layout = "vertical",
              verticalAlign="middle") %>% 
    hc_colorAxis(  stops= cor_colr,min=rate,max=1) %>%
    dapar_hc_ExportMenu(filename = "corrMatrix")
}