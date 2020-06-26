
#' @title Builds a plot from a numeric matrix.
#' 
#' @param qDataBefore A numeric matrix that contains quantitative data before normalization.
#' 
#' @param qDataAfter A numeric matrix that contains quantitative data after normalization.
#' 
#' @param conds A vector of the conditions (one condition per sample).
#' 
#' @param palette xxx
#' 
#' @param subset.view xxx
#' 
#' @param n xxx
#' 
#' @param scatter xxx
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(highcharter)
#' library(Features)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' conds <- colData(obj)[["Condition"]]
#' obj <- normalizeD(obj, 2, name='norm', method='SumByColumns', conds=conds, type='overall')
#' compareNormalizationD_HC(assay(obj, 2), assay(obj, 3), conds=conds, n=100)
#' 
#' @import highcharter
#' 
#' @export
#' 
compareNormalizationD_HC <- function(qDataBefore,
                                     qDataAfter,
                                     conds =NULL,
                                     palette = NULL,
                                     subset.view = NULL,
                                     n = NULL,
                                     type = 'scatter'){
  
  if (is.null(conds)){
    warning("'conds' is null.")
    return(NULL)
  }
  
  print(str(subset.view))
  print(paste0('subset.view:',subset.view))
  
  if (!is.null(subset.view) && length(subset.view) > 0)
  {
    if (nrow(qDataBefore) > 1)
      if (length(subset.view)==1){
        qDataBefore <- as_tibble(cbind(t(qDataBefore[subset.view,])))
        qDataAfter <- as_tibble(cbind(t(qDataAfter[subset.view,])))
      } else {
        qDataBefore <- as_tibble(cbind(qDataBefore[subset.view,]))
        qDataAfter <- as_tibble(cbind(qDataBefore[subset.view,]))
      }
  }
  
  
  if (!match(type, c('scatter', 'line') )){
    warning("'type' must be equal to 'scatter' or 'line'.")
    return(NULL)
  }
  
  
  if (is.null(n)){
    n <- seq_len(nrow(qDataBefore))
  } else {
    if (n > nrow(qDataBefore)){
      warning("'n' is higher than the number of rows of datasets. Set to number 
              of rows.")
      n <- nrow(qDataBefore)
    }
    
    ind <- sample(seq_len(nrow(qDataBefore)),n)
    if (nrow(qDataBefore) > 1)
      if (length(ind) == 1){
        qDataBefore <- as_tibble(cbind(t(qDataBefore[ind,])))
        qDataAfter <- as_tibble(cbind(t(qDataAfter[ind,])))
      } else {
        qDataBefore <- as_tibble(cbind(qDataBefore[ind,]))
        qDataAfter <- as_tibble(cbind(qDataAfter[ind,]))
      }
  }
  
  palette <- BuildPalette(conds, palette)
  

  x <- qDataBefore
  y <- qDataAfter/qDataBefore
  
  ##Colors definition
  legendColor <- unique(palette)
  txtLegend <- unique(conds)
  
  
  series <- list()
  for (i in 1:length(conds)){
    tmp <- list(name=colnames(x)[i], data =list_parse(data.frame(x=x[,i],
                                                                 y=y[,i])
                                                      )
    )
    series[[i]] <- tmp
  }
  
  h1 <-  highchart() %>% 
    dapar_hc_chart( chartType = type) %>%
    hc_add_series_list(series) %>%
    hc_colors(palette) %>%
    hc_tooltip(enabled= "false" ) %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  h1
  
}

