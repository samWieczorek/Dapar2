
#' @title Builds a complete color palette for the conditions given in argument
#' 
#' @description xxxx
#' 
#' @param conds The extended vector of samples conditions
#' 
#' @param base_palette The basic color (HEX code) used to build the complete palette. This vector have the same length as unique(conds)
#' 
#' @return A vector composer of HEX color code for the conditions
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' conditions <- SummarizedExperiment::colData(Exp1_R25_pept)$Condition
#' BuildPalette(conditions, c('AAAAAA', 'BBBBBB'))
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' 
BuildPalette <- function(conds, base_palette){
  
  palette <- NULL
  if (is.null(base_palette)){
    palette.init <- RColorBrewer::brewer.pal(8,"Dark2")[1:max(3,length(unique(conds)))]
  } else {
    palette.init <- base_palette
  }
  
  for (i in 1:length(conds)){
    palette[i] <- palette.init[which(conds[i] == unique(conds))]
  }
  return(palette)
  
}                                        




#' @title Returns the number of empty lines in the data
#' 
#' @param qData A matrix corresponding to the quantitative data.
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept,1)
#' nEmptyLines(qData)
#' 
#' @export
#' 
nEmptyLines <- function(qData){
  n <- sum(apply(is.na(as.matrix(qData)), 1, all))
  return (n)
}





#' @title Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#' 
#' @param data A data.frame
#' 
#' @param type The value to search in the dataframe
#' 
#' @return A boolean dataframe 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' data <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$OriginOfValues]
#' is.OfType(as.data.frame(data@listData), "MEC")
#' 
#' @export
#' 
is.OfType <- function(data, type){
  
  if (class(data) != 'data.frame'){
    warning('The parameter param is not a data.frame.')
    return(NULL)
  }
  return (type == data)
}






#' @title Similar to the function \code{is.na} but focused on the equality with the missing 
#' values in the dataset (type 'POV' and 'MEC')
#' 
#' @param data A data.frame
#' 
#' @return A boolean dataframe 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' data <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$OriginOfValues]
#' is.MV(as.data.frame(data@listData))
#' 
#' @export
#' 
is.MV <- function(data){
  if (class(data) != 'data.frame'){
    warning('The parameter param is not a data.frame.')
    return(NULL)
  }
  POV = is.OfType(data, "POV")
  MEC = is.OfType(data, "MEC")
  isNA = is.na(data)
  df <- POV | MEC | isNA
  
  return (df)
}




#' @title Returns the possible number of values in lines in the data
#' 
#' @param obj An object of class \code{QFeatures}
#' 
#' @param i The indice of the dataset (SummarizedExperiment) in the object
#' 
#' @param type wholeMatrix, allCond or atLeastOneCond
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' getListNbValuesInLines(Exp1_R25_pept, 1)
#' 
#' @export
#' 
#' @importFrom SummarizedExperiment colData rowData 
#' 
#' @importFrom S4Vectors sort
#' 
getListNbValuesInLines <- function(obj, i, type="wholeMatrix"){
  
  if (is.null(obj)){return(NULL)}
  
  if(is.null(metadata(obj)$OriginOfValues)){
    data <- as.data.frame(assay(obj[[i]]))
  }
  else {
    data <- as.data.frame(rowData(obj[[i]])[, metadata(obj)$OriginOfValues])
  }
  
  switch(type,
         wholeMatrix= {
           ll <- unique(ncol(data) - apply(is.MV(data), 1, sum))
         },
         allCond = {
           tmp <- NULL
           for (cond in unique(colData(obj)[['Condition']])){
             tmp <- c(tmp, length(which(colData(obj)[['Condition']] == cond)))
           }
           ll <- seq(0,min(tmp))
         },
         atLeastOneCond = {
           tmp <- NULL
           for (cond in unique(colData(obj)[['Condition']])){
             tmp <- c(tmp, length(which(colData(obj)[['Condition']] == cond)))
           }
           ll <- seq(0,max(tmp))
         }
  )
  
  return (sort(ll))
}





#' @title Customised contextual menu of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param filename The filename under which the plot has to be saved
#' 
#' @return A contextual menu for highcharts plots
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc_chart(hc,type = "line") 
#' hc_add_series(hc,data = c(29, 71, 40))
#' dapar_hc_ExportMenu(hc,filename='foo')
#' 
#' @export
#' 
#' @importFrom highcharter hc_exporting
#' 
dapar_hc_ExportMenu <- function(hc, filename){
  hc_exporting(hc, enabled=TRUE,
               filename = filename,
               buttons= list(
                 contextButton= list(
                   menuItems= list('downloadPNG', 'downloadSVG','downloadPDF')
                 )
               )
  )
}




#' @title Customised resetZoomButton of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param chartType The type of the plot
#' 
#' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
#' 
#' @param width xxx
#' 
#' @param height xxx
#' 
#' @return A highchart plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc_chart(hc,type = "line") 
#' hc_add_series(hc,data = c(29, 71, 40))
#' dapar_hc_ExportMenu(hc,filename='foo')
#' 
#' @export
#' 
#' @importFrom highcharter hc_chart
#' 
dapar_hc_chart <- function(hc,  chartType, zoomType="None", width=0, height=0){
  hc %>% 
    hc_chart(type = chartType, 
             zoomType=zoomType,
             showAxes = TRUE,
             width = width,
             height = height,
             resetZoomButton= list(
               position = list(
                 align= 'left',
                 verticalAlign = 'top')
             ))
}




#' @title Retrieve the indices of non-zero elements in sparse matrices
#' 
#' @description This function retrieves the indices of non-zero elements in sparse matrices
#' of class dgCMatrix from package Matrix. This function is largely inspired from 
#' the package \code{RINGO}
#' 
#' @param x A sparse matrix of class dgCMatrix
#' 
#' @return A two-column matrix
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(Matrix)
#' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, sparse=TRUE)
#' res <- nonzero(mat)
#' 
#' @export
#' 
nonzero <- function(x){
  ## function to get a two-column matrix containing the indices of the
  ### non-zero elements in a "dgCMatrix" class matrix
  
  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0))
    return(matrix(0, nrow=0, ncol=2,
                  dimnames=list(character(0), c("row","col"))))
  res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  colnames(res) <- c("row", "col")
  res <- res[x@x != 0, , drop = FALSE]
  return(res)
}

