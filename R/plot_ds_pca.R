
#' @param object An instance of the class `SummarizedExperiment`
#' @param conds xxx
#' @param var.scaling The dimensions to plot
#' @param ncp A `integer(1)` which represents the umber of dimensions kept in the results.
#' 
#' @importFrom FactoMineR PCA
#' @importFrom stats na.omit
#' 
#' @export
#' 
#' @rdname descriptive-statistics
#' 
wrapper_pca <- function(object, 
                        conds, 
                        var.scaling = TRUE, 
                        ncp = NULL){
  
  if(missing(object))
    stop("'object' is missing.")
  if(missing(conds))
    stop("'conds' is missing.")
  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  
  qData <- assay(object)
  
  if (is.null(var.scaling))
    var.scaling <- TRUE
  
  if (length(which(is.na(qData))) > 0)
    qData <- stats::na.omit(qData)
  
  
  if (is.null(ncp)){
    nmax <- 12
    y <- qData
    nprot <- dim(y)[1]
    n <- dim(y)[2] # If too big, take the number of conditions.
    
    if (n > nmax){
      n <- length(unique(conds))
    }
    
    ncp <- min(n, nmax)
  }
  
  res.pca <- FactoMineR::PCA(qData, 
                 scale.unit = var.scaling, 
                 ncp = ncp, 
                 graph= FALSE
                 )
  
  return(res.pca)
}



#' #' @param res.pca Result of FactoMineR::PCA
#' #' @param chosen.axes The dimensions to plot
#' #' 
#' #' @author Samuel Wieczorek, Enora Fremy
#' #' 
#' #' @importFrom factoextra fviz_pca_var
#' #' 
#' #' @export
#' #' 
#' #' @rdname pca_plots
#' #' 
#' plotPCA_Var <- function(res.pca = NULL, 
#'                         chosen.axes = c(1, 2)){
#'   if (is.null(res.pca))
#'     return(NULL)
#'   
#'   fviz_pca_var(res.pca, 
#'                axes = chosen.axes,
#'                col.var = "cos2",
#'                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#'                repel = TRUE
#'                )
#' }



#' #' @param res.pca Result of FactoMineR::PCA
#' #' @param axes The dimensions to plot
#' #' @param geom xxx
#' #' 
#' #' @author Samuel Wieczorek, Enora Fremy
#' #' 
#' #' @importFrom factoextra fviz_pca_ind
#' #' 
#' #' @export
#' #' 
#' #' @rdname pca_plots
#' #' 
#' plotPCA_Ind <- function(res.pca, chosen.axes = c(1,2)){
#'   if (is.null(res.pca))
#'     return(NULL)
#'   
#'   fviz_pca_ind(res.pca,
#'                axes = chosen.axes,
#'                geom  ="point")
#'   
#' }





#' @param res.pca Result of FactoMineR::PCA
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @export
#' 
#' @import highcharter
#' 
#' @rdname pca_plots
#' 
plotPCA_Eigen <- function(res.pca){
  stopifnot(!is.null(res.pca))
  
  hc <- highchart() %>%
    hc_yAxis_multiples(list(title = list(text = "% of variances"),
                            lineWidth = 0, 
                            labels = list(format = "{value}%"), 
                            max = 100), 
                       list(title = list(text = "Cumulative % of variances"), 
                            opposite = FALSE, 
                            max = 100),
                       list(title = list(text = "Eigen values"), 
                            opposite = TRUE, 
                            labels = list(format = "{value}%") )
                       ) %>%
    hc_xAxis(title = "Principal Components", categories = rownames(res.pca$eig)) %>%
    hc_add_series(data.frame(y = res.pca$eig[,2]), 
                  type = "column",  
                  name = "% of variances", 
                  yAxis = 0) %>%
    hc_add_series(data.frame(y = res.pca$eig[,3]), 
                  type = "line",  
                  color="darkblue",
                  name = "Cumulative % of variances", 
                  marker = "diamond", 
                  color = "#FF7900", 
                  yAxis = 0) %>%
    hc_legend(enabled = TRUE)
  
  hc
}
