#' @export 
#' @import highcharter
#' @rdname descriptive-statistics
mvPerLinesHisto <- function(object){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  qData <- assay(object)
  
  coeffMax <- .1
  
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  
  nb.col <- dim(qData)[2] 
  nb.na <- NbNAPerRow
  temp <- table(NbNAPerRow)
  nb.na2barplot <- c(temp, rep(0,1+ncol(qData)-length(temp)))
  
  if (sum(NbNAPerRow) == 0){
    nb.na2barplot <- rep(0, 1 + ncol(qData))
  }
  
  df <- data.frame(y=nb.na2barplot[-1])
  myColors = rep("lightgrey",nrow(df))
  myColors[nrow(df)] <- "red"
  
  
  
  h1 <-  highchart() %>% 
    hc_title(text = "#[lines] with X NA values") %>% 
    hc_add_series(data = df, type="column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation = list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(df), title = list(text = "#[NA values] per line")) %>%
    customExportMenu(fname = "mvPlot1") %>%
    hc_tooltip(enabled = TRUE,
               headerFormat= '',
               pointFormat = "{point.y} ")
  
  return(h1)
  
}



#' @export
#' 
#' @import highcharter
#' @rdname descriptive-statistics
mvPerLinesHistoPerCondition <- function(object, 
                                        conds, 
                                        pal.name){
  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  qData <- assay(object)
  
  myColors <-  SampleColors(conds, pal.name)
  
  nbConditions <- length(unique(conds))
  
  ncolMatrix <- max(unlist(lapply(unique(conds), function(x){length(which(conds==x))})))
  m <- matrix(rep(0, nbConditions*(1+ncolMatrix)), 
              ncol = nbConditions, 
              dimnames=list(seq(0:(ncolMatrix)), unique(conds)))
  
  for (i in unique(conds)) {
    nSample <- length(which(conds == i))
    t <- NULL
    if (nSample == 1) {
      t <- table(as.integer(is.na(qData[,which(conds == i)])))
    } else {t <- table(rowSums(is.na(qData[ ,which(conds == i)])))}
    
    m[as.integer(names(t))+1,i] <- t
  }
  
  m <- as.data.frame(m)
  
  rownames(m) <- 0:(nrow(m)-1)
  
  h1 <-  highchart() %>% 
    hc_title(text = "#[lines] with X NA values (condition-wise)") %>% 
    customChart(chartType = "column") %>%
    hc_plotOptions( column = list(stacking = ""),
                    dataLabels = list(enabled = FALSE),
                    animation = list(duration = 100)) %>%
    hc_colors(unique(myColors)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(m), title = list(text = "#[NA values] per line (condition-wise)")) %>%
    customExportMenu(fname = "mvPlot_2") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y} ")
  
  for (i in seq_len(nbConditions)){
    h1 <- h1 %>% hc_add_series(data=m[ ,unique(conds)[i]]) }
  
  
  return(h1)
  
}




#' @export
#' @import highcharter
#' @rdname descriptive-statistics
mvHisto <- function(object, 
                       conds, 
                       showValues = FALSE, 
                       pal.name = NULL){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  qData <- assay(object)
  
  myColors <-  SampleColors(conds, pal.name)
  NbNAPerCol <- colSums(is.na(qData))
  NbNAPerRow <- rowSums(is.na(qData))
  
  df <- data.frame(NbNAPerCol)
  names(df) <- 'y'
  
  
  h1 <-  highchart() %>%
    customChart(chartType = "column") %>%
    hc_title(text = "#NA by replicate") %>%
    hc_add_series(df,type="column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds, title = list(text = "Replicates")) %>%
    customExportMenu(fname = "mvPlot_3") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y}")
  
  
  return(h1)
  
}


