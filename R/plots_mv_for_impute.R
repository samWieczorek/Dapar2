#' @title Distribution of Observed values with respect to intensity values
#' 
#' @description 
#' 
#' This method shows density plots which represents the repartition of
#' Partial Observed Values for each replicate in the dataset.
#' The colors correspond to the different conditions (Condition in colData of
#' the object of class [QFeatures]).
#' The x-axis represent the mean of intensity for one condition and one
#' entity in the dataset (i. e. a protein) 
#' whereas the y-axis count the number of observed values for this entity
#' and the considered condition.
#' 
#' @section xxx:
#' 
#' xxxx
#' 
#' @section missing value image:
#' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class \code{MSnSet} and 
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to visualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' 
#' @name imputation

#' @return A plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' 
#' library(QFeatures)
#' data(ft)
#' qData <- assay(ft,1)
#' conds <- design(ft)$Condition
#' 
#' #-----------------------------
#' # xxx
#' #-----------------------------
#' 
#' mvPlot.2(qData, conds)
#' 
#' mvPlot.2(Exp1_R25_pept, 2, title="POV distribution", palette=pal)
#' 
#' #-----------------------------
#' # xxx
#' #-----------------------------
#' 
#' data(ft_na)
#' mvImage(ft_na[[1]], conds)
#' 
#' #-----------------------------
#' # xxx
#' #-----------------------------
#' 
#' if(interactive()){
#'  library(QFeatures)
#'  library(shiny)
#'  library(DaparToolshed)
#'  data(ft_na)
#'  conds <- design(ft_na)$Condition
#'  
#'   
#'  ui <- mod_imputation_plot_ui('plot')
#' 
#'  server <- function(input, output, session) {
#'   mod_imputation_plot_server('plot',
#'                           obj = reactive({ft_na[[1]]}),
#'                           conds = reactive({conds}),
#'                           pal.name = reactive({'Dark2'})
#'                          )
#'   }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
NULL


#' @param qData An instance of a `matrix` containing numeric vales.
#' @param conds A `character()` of the name of conditions (one condition per sample).
#' @param title The title of the plot
#' @param pal.name A `character(1)` which is the name of the palette (from
#' the package [RColorBrewer] to use.
#' 
#' @export
#' 
#' @import highcharter
#' 
#' @rdname imputation
#' 
mvPlot.2 <- function(qData,
                     conds,
                     title = NULL,
                     pal.name){
  
  if(missing(qData))
    stop("'qData' is required.")
  
  stopifnot(inherits(qData, 'matrix'))
  
  if(missing(conds))
    stop("'conds' is required.")
  
  if(missing(pal.name))
    pal.name <- NULL
  
  myColors <- NULL
  uconds <- unique(conds)
  palette <- SampleColors(uconds, pal.name)
  mTemp <- nbNA <- nbValues <- matrix(rep(0,nrow(qData)*length(uconds)), 
                                      nrow = nrow(qData),
                                      dimnames = list(NULL, uconds)
                                      )
  dataCond <- data.frame()
  ymax <- 0
  series <- list()
  myColors <- NULL
  j <- 1
  
  for (iCond in uconds){
    if (length(which(conds==iCond)) == 1){
      
      mTemp[,iCond] <- qData[,which(conds==iCond)]
      nbNA[,iCond] <- as.integer(is.OfType(qData[,which(conds==iCond)]))
      nbValues[,iCond] <- length(which(conds==iCond)) - nbNA[,iCond]
    } else {
      mTemp[,iCond] <- apply(qData[,which(conds==iCond)], 1, mean, na.rm=TRUE)
      nbNA[,iCond] <- apply(qData[,which(conds==iCond)],1,function(x) length(which(is.na(x) == TRUE)))
      nbValues[,iCond] <- length(which(conds==iCond)) - nbNA[,iCond]
    }
    
    
    for (i in 1:length(which(conds==iCond))){
      data <- mTemp[which(nbValues[, iCond] == i), iCond]
      tmp <- NULL    
      if (length(data) >= 2)
      {
        tmp <- density(mTemp[which(nbValues[,iCond]==i),iCond])
        tmp$y <- tmp$y + i
        if (max(tmp$y) > ymax) { ymax <- max(tmp$y)}
      }
      series[[j]] <- tmp
      myColors <- c(myColors, palette[which(unique(conds)==iCond)])
      j <- j+1
    }
    
  }
  
  axis.title <- list(x = "Mean of intensities",
                     y = "Number of quantity values per condition")
  
  hc <-  highchart() %>%
    hc_title(text = title) %>%
    customChart(chartType = "spline", zoomType="xy") %>%
    
    hc_legend(align = "left", 
              verticalAlign = "top",
              layout = "vertical") %>%
    hc_xAxis(title = list(text = axis.title$x)) %>%
    hc_yAxis(title = list(text = axis.title$y),
             tickInterval= 1
             ) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    customExportMenu(fname = "POV_distribution") %>%
    hc_plotOptions(
      series=list(
        showInLegend = TRUE,
        connectNulls= TRUE,
        marker=list(
          enabled = FALSE)
        
      )
    )
  
  for (i in 1:length(series)){
    hc <- hc_add_series(hc,
                        data = list_parse(data.frame(cbind(x = series[[i]]$x, 
                                                           y = series[[i]]$y))), 
                        showInLegend=FALSE,
                        color = myColors[i],
                        name=conds[i])
  }
  
  # add three empty series for the legend entries. Change color and marker symbol
  for (c in 1:length(uconds)){
    hc <-  hc_add_series(hc,data = data.frame(),
                         name = uconds[c],
                         color = palette[c],
                         marker = list(symbol = "circle"),
                         type = "line")
  }
  
  return(hc)
}



#' @param object An instance of a class `SummarizedExperiment`.
#' @param conds A `character()` of the name of conditions (one condition per sample).
#' 
#' @rdname imputation
#' @export
#' @importFrom stats setNames
mvImage <- function(object,
                    conds
                    ){
  
  if(missing(object))
    stop("'object' is required.")
  stopifnot(inherits(object, "SummarizedExperiment"))
  
  if(missing(conds))
    stop("'conds' is required.")

  stopifnot('qMetadata' %in% colnames(rowData(object)))
  
  i_MEC <- which(apply(qMetadata(object)=="missing MEC", 1, sum) >0)
  
  if (length(i_MEC)==0){
    warning("The dataset contains no Missing value on Entire Condition (MEC). 
            So this plot is not available.")
    return(NULL)
  }else if (length(i_MEC)==1){
    warning("The dataset contains only one Missing value on Entire Condition (MEC). 
            Currently, Prostar does not handle such dataset to build the plot. 
            As it has no side-effects on the results, you can continue your imputation.")
    return(NULL)
  }
  
  qData <- assay(object)[i_MEC,]
  
  
  ### build indices of conditions
  indCond <- list()
  conds.names <- unique(conds)
  
  indCond <- setNames(lapply(conds.names, 
                             function(x) 
                               grep(x, conds)), 
                      paste0('cond', seq_len(length(conds.names))))
  
  nNA1 = apply(as.matrix(qData[,indCond$cond1]), 
               1, 
               function(x) 
                 sum(is.na(x))
               )
  nNA2 = apply(as.matrix(qData[,indCond$cond2]), 
               1, 
               function(x) 
                 sum(is.na(x))
               )
  
  o <- order(((nNA1 + 1)^2) / (nNA2 + 1))
  exprso <- qData[o,]
  
  for (i in 1:nrow(exprso)){
    k <- order(exprso[i,indCond$cond1])
    exprso[i, rev(indCond$cond1)] <- exprso[i, k]
    .temp <- mean(exprso[i, rev(indCond$cond1)], na.rm = TRUE)
    exprso[i, which(!is.na(exprso[i,indCond$cond1]))] <- .temp
    
    k <- order(exprso[i,indCond$cond2])
    exprso[i, indCond$cond2] <- exprso[i, k+length(indCond$cond1)]
    .temp <- mean(exprso[i, indCond$cond2], na.rm = TRUE)
    exprso[i, length(indCond$cond1) + 
             which(!is.na(exprso[i, indCond$cond2]))] <- .temp
  }
  
  colors <- colorRampPalette(c("yellow", "red"))(100)
  mv.heatmap(exprso,
                col = colors,
                key=TRUE,
                srtCol= 0,
                labCol=conds,
                ylab = "Peptides / proteins",
                main = "MEC heatmap"
  )
  
  #heatmap_HC(exprso,col = colfunc(100),labCol=conds)
}

