#' Boxplot for quantitative proteomics data using the library \code{highcharter}
#' 
#' @title Builds a boxplot from a dataframe using the library \code{highcharter}
#' @param qData Numeric matrix 
#' @param fData xxx
#' @param legend A vector of the conditions (one condition per sample).
#' @param palette xxx
#' @param subset.view A vector of index indicating rows to highlight
#' @return A boxplot
#' @author Samuel Wieczorek, Anais Courtier, Enora Fremy
#' @seealso \code{\link{densityPlotD_HC}}
#' @examples
#' library(Features)
#' library(DAPAR2)
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original_log']])
#' fData <- rowData(Exp1_R25_pept[['original_log']])[['Sequence']]
#' legend <- colData(Exp1_R25_pept)[["Sample.name"]]
#' boxPlotD_HC(qData, fData, legend, subset.view=1:10)
#' @importFrom Biobase exprs pData
#' @import highcharter
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
boxPlotD_HC <- function(qData, fData, legend=NULL, palette = NULL,subset.view=NULL){
  
  if (is.null(qData)){
    warning('The dataset in NULL and cannot be shown')
    return(NULL)
  }
  
  if( is.null(legend)){legend <- colData(Exp1_R25_pept)[["Sample.name"]] }
  
  palette <- BuildPalette(colData(Exp1_R25_pept)[["Condition"]], palette)
  
  bx <-boxplot(qData, na.rm=TRUE)
  df_outlier <- data.frame(x=bx$group-1,
                           y = bx$out)
  
  #tmp <- NULL
  #for (i in 1:ncol(qData)){
  #  tmp <- c(tmp, rep(paste("S",i,legend[i], collapse="_"),nrow(qData)))
  #}
  tmp <- NULL
  for (i in 1:ncol(qData)){
    tmp <- c(tmp, rep(paste(paste0(rep("A", i), collapse=""),legend[i], sep='_'),nrow(qData)))
  }
  
  df <- data.frame(values = as.vector(qData,mode='numeric'),samples = tmp, stringsAsFactors = FALSE)
  
  
  hc <- highcharter::hcboxplot(x=df$values, var = df$samples, colorByPoint = TRUE, outliers = TRUE) %>%
    highcharter::hc_chart(type="column") %>%
    highcharter::hc_yAxis(title = list(text = "Log (intensity)")) %>%
    highcharter::hc_xAxis(title = list(text = "Samples"), categories=legend) %>%
    highcharter::hc_colors(palette) %>%
    highcharter::hc_add_series(type= "scatter",df_outlier,name="Outliers",tooltip=list(enabled=F,headerFormat ="",pointFormat="{point.y: .2f} ")) %>%  
    highcharter::hc_plotOptions(
      
      boxplot= list(
        
        fillColor= c('lightgrey'),
        lineWidth= 3,
        medianColor= 'grey',
        medianWidth= 3,
        stemColor= '#A63400',
        stemDashStyle= 'dot',
        stemWidth= 1,
        whiskerColor= '#3D9200',
        whiskerLength= '20%',
        whiskerWidth= 3
      ),
      scatter = list(
        marker=list(
          fillColor = 'white',
          lineWidth = 0.5,
          lineColor = 'grey'
        )
      )
    )
  
  # Display of rows to highlight (index of row in subset.view) 
  if(!is.null(subset.view)){
    idVector <- fData
    pal=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(subset.view))    
    n=0
    for(i in subset.view){
      n=n+1
      dfSubset <- data.frame(y = as.vector(qData[i,],mode='numeric'), x = as.numeric(factor(names(qData[i,])))-1, stringsAsFactors = FALSE)
      hc<-hc %>%
        highcharter::hc_add_series(type= "line",data=dfSubset,color=pal[n], dashStyle = "shortdot",name=idVector[i],
                                   tooltip=list(enabled=T,headerFormat ="",pointFormat="{point.series.name} : {point.y: .2f} ") )
      
    }
  }  
  
  hc
  
}