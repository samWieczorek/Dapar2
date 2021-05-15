
#' @title Returns the number of empty lines in the data
#' 
#' @param qData A matrix corresponding to the quantitative data.
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
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
#' library(QFeatures)
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
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' data <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$OriginOfValues]
#' is.MV(as.data.frame(data@listData))
#' 
#' @export
#' 
is.MV <- function(data){
  # if (class(data) != 'data.frame'){
  #   warning('The parameter param is not a data.frame.')
  #   return(NULL)
  # }
  #POV = is.OfType(data, "POV")
  #MEC = is.OfType(data, "MEC")
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
#' @param type WholeMatrix, AllCond or AtLeastOneCond
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
getListNbValuesInLines <- function(obj, i, type="WholeMatrix"){
  
  if (is.null(obj)){return(NULL)}
  if(missing(i))
    stop("'i' is required")
  if (!(type %in% c('None', 'WholeMatrix', 'AtLeastOneCond', 'AllCond')))
    stop("'type' is not one of: 'None', 'WholeMatrix', 'AtLeastOneCond', 'AllCond'")
  
  data <- as.data.frame(assay(obj[[i]]))
  
  switch(type,
         WholeMatrix= {
           ll <- unique(ncol(data) - apply(is.na(data), 1, sum))
         },
         AllCond = {
           tmp <- NULL
           for (cond in unique(colData(obj)[['Condition']])){
             tmp <- c(tmp, length(which(colData(obj)[['Condition']] == cond)))
           }
           ll <- seq(0,min(tmp))
         },
         AtLeastOneCond = {
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

