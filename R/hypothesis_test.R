

#' @title Record the adjacency matrices 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' Exp1_R25_pept <- addAssay(Exp1_R25_pept, QFeatures::filterNA(Exp1_R25_pept[[2]],  pNA = 0), name='filtered')
#' 
#' Exp1_R25_pept <- HypothesisTest(object = Exp1_R25_pept, i = 3, test = 't_test', contrast='OnevsOne', type='Student')
#' 
#' Exp1_R25_pept <- HypothesisTest(object = Exp1_R25_pept, i=3, test = 'limma', comp.type = "OnevsOne")
#' 
#' Exp1_R25_pept <- HypothesisTest(object = Exp1_R25_pept, i=3, test = '1wayAnova', with_post_hoc = TRUE, post_hoc_test = 'Dunnett')
#' 
#' @export
#' 
"HypothesisTest"


#' @param  obj.se An object (peptides) of class \code{SummarizedExperiment}.
#' 
#' @return The slot processing of obj@processingData
#' 
setMethod('HypothesisTest', "SummarizedExperiment",
          function(object, test, sTab, ...) {
            
            res.ll <- switch(test,
                             t_test = compute_t_tests(obj = object, 
                                                      sTab = sTab,
                                                      ...),
                             limma = limmaCompleteTest(obj = object,
                                                       sTab = sTab,
                                                       ...),
                             '1wayAnova' = wrapperClassic1wayAnova(obj = object,
                                                                   sTab = sTab,
                                                                   ...)
            )
            
            rowData(object) <-  cbind(rowData(object), res.ll)
            
            metadata(object)$HypothesisTest_cols <- colnames(res.ll)
            object
          }
)


setMethod("HypothesisTest", "QFeatures",
          function(object, 
                   i, 
                   name = "HypothesisTestAssay", 
                   test = NULL,
                   ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            
            
            object <- addAssay(object,
                               HypothesisTest(object[[i]], 
                                              test = test,
                                              sTab = colData(object),
                                              ...
                                              ),
                               name)
            addAssayLink(object, 
                         from = i, 
                         to = name)
            
            
          }
)

























#' @title Density plots of logFC values
#' 
#' @description This function show the density plots of Fold Change (the same as calculated by limma) for a list 
#' of the comparisons of conditions in a differential analysis.
#' 
#' @param df_logFC A dataframe that contains the logFC values for the set of comparisons to show. Each column corresponds to a comparison.
#' 
#' @param threshold_LogFC The threshold on log(Fold Change) to distinguish between differential and non-differential data 
#' 
#' @param palette xxx
#' 
#' @return A highcharts density plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' obj <- addAssay(obj, QFeatures::filterNA(obj[[2]],  pNA = 0), name='filtered')
#' se <- t_test_sam(obj[[3]], colData(obj), FUN = compute.t.test)
#' ind <- grep('_logFC', colnames(metadata(se)$t_test))
#' df <- setNames(as.data.frame(metadata(se)$t_test[,ind]), colnames(metadata(se)$t_test)[ind])
#' pal <- ExtendPalette(2, 'Dark2')
#' if (length(ind)>0)
#' hc_logFC_DensityPlot(df, palette = pal)
#' 
#' @export
#' 
#' @import highcharter
#' 
hc_logFC_DensityPlot <-function(df_logFC, 
                                threshold_LogFC = 0, 
                                palette=NULL){
  
  if (is.null(df_logFC) || threshold_LogFC < 0){
    hc <- NULL
    return(NULL)
  }
  
  
  myColors <- NULL
  if (is.null(palette)){
    warning("Color palette set to default.")
    palette <- ExtendPalette(ncol(df_logFC), "Paired")
  } else {
    if (length(palette) != ncol(df_logFC)){
      warning("The color palette has not the same dimension as the number of samples")
      palette <- ExtendPalette(ncol(df_logFC), "Paired")
    }
  }
  
  nValues <- nrow(df_logFC)*ncol(df_logFC)
  nInf <- length(which(df_logFC <= -threshold_LogFC))
  nSup <- length(which(df_logFC >= threshold_LogFC))
  nInside <- length(which(abs(df_logFC) < threshold_LogFC))
  
  hc <-  highcharter::highchart() %>% 
    hc_title(text = "log(FC) repartition") %>% 
    dapar_hc_chart(chartType = "spline", zoomType="x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_colors(palette) %>%
    hc_xAxis(title = list(text = "log(FC)"),
             plotBands = list(list(from= -threshold_LogFC, to = threshold_LogFC, color = "lightgrey")),
             plotLines=list(list(color= "grey" , width = 2, value = 0, zIndex = 5)))%>%
    hc_yAxis(title = list(text="Density")) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "densityplot") %>%
    hc_plotOptions(
      series=list(
        animation=list(duration = 100),
        connectNulls= TRUE,
        marker=list(enabled = FALSE)
      )
    )
  
  maxY.inf <- NULL
  maxY.inside <- NULL
  maxY.sup <- NULL
  minX <- NULL
  maxX<- NULL
  
  
  for (i in 1:ncol(df_logFC)){
    tmp <- density(df_logFC[,i],na.rm = TRUE)
    ind <- tmp$y[which(tmp$x <= -threshold_LogFC)]
    maxY.inf <- max(maxY.inf, ifelse(length(ind)==0,0,ind))
    maxY.inside <- max(maxY.inf, tmp$y[intersect(which(tmp$x > -threshold_LogFC),which(tmp$x < threshold_LogFC))])
    ind <- tmp$y[which(tmp$x > threshold_LogFC)]
    maxY.sup <- max(maxY.sup, ifelse(length(ind)==0,tmp$y[length(tmp$y)],ind))
    minX <- min(minX,tmp$x)
    maxX <- max(maxX,tmp$x)
    
    
    hc <- hc_add_series(hc,
                        data.frame(x = tmp$x,  y = tmp$y), 
                        name=colnames(df_logFC)[i])
  }
  
  
  if(threshold_LogFC != 0) {
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape='connector',
        backgroundColor = 'lightgrey',
        #verticalAlign = 'bottom',
        align='left',
        #distance=0,
        style=list(
          fontSize= '1.5em',
          textOutline= '1px white'
        ),
        borderWidth = 0,
        x = 20
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = 0,
            y = maxY.inside
          ),
          text = paste0("n Filtered out = ",nInside, "<br>(", round(100*nInside/nValues, digits=2), "%)")
        )
      )
    )
  }
  if (threshold_LogFC >= minX){
    hc <- hc %>%
      hc_add_annotation(
        labelOptions = list(
          shape='connector',
          backgroundColor = 'rgba(255,255,255,0.5)',
          verticalAlign = 'top',
          borderWidth = 0,
          crop=TRUE,
          style=list(
            color = 'blue',
            fontSize= '1.5em',
            textOutline= '1px white'
          ),
          y = -10
        ),
        labels = list(
          list(
            point = list(
              xAxis = 0,
              yAxis = 0,
              x = mean(c(minX,-threshold_LogFC)),
              y = maxY.inf
            ),
            text = paste0("nInf = ",nInf, "<br>(", round(100*nInf/nValues, digits=2), ")%")
          )
        )
      )
  }
  
  if (threshold_LogFC <= maxX){
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape='connector',
        backgroundColor = 'blue',
        verticalAlign = 'top',
        borderWidth = 0,
        style=list(
          color = 'blue',
          fontSize= '1.5em',
          textOutline= '1px white'
        ),
        y = -5
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x =  mean(c(maxX,threshold_LogFC)),
            y = maxY.sup
          ),
          text = paste0("nSup = ",nSup, "<br>(", round(100*nSup/nValues, digits=2), ")%")
        )
      )
    )
    
  }
  
  return(hc)
  
}

