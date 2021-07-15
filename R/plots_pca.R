


#' @title Compute the PCA
#' 
#' @param qData numeric matrix
#' 
#' @param condition xxx
#' 
#' @param var.scaling The dimensions to plot
#' 
#' @param ncp Number of dimensions kept in the results
#' 
#' @return A list including eigenvalues of obj
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' condition <- colData(Exp1_R25_pept)[["Condition"]]
#' res.pca <- wrapper.pca(qData, condition)
#' 
#' plotPCA_Var(res.pca)
#' 
#' plotPCA_Ind(res.pca)
#' 
#' plotPCA_Eigen_hc(res.pca)
#' 
#' 
"wrapper.pca"

#' @importFrom FactoMineR PCA
#' @importFrom stats na.omit
#' 
#' @export
#' 
wrapper.pca <- function(qData, 
                        condition, 
                        var.scaling = TRUE, 
                        ncp = NULL){
  
  if(missing(qData))
    stop("'qData' is missing.")
  if(missing(condition))
    stop("'condition' is missing.")
  
  if (is.null(var.scaling)){
    warning("Setting 'var.scaling' to TRUE.")
    var.scaling <- TRUE
  }
  
  if (length(which(is.na(qData))) > 0)
    qData <- stats::na.omit(qData)
  
  
  if (is.null(ncp)){
    nmax <- 12
    y <- qData
    nprot <- dim(y)[1]
    n <- dim(y)[2] # If too big, take the number of conditions.
    
    if (n > nmax){
      n <- length(unique(condition))
    }
    
    ncp <- min(n, nmax)
  }
  
  res.pca <- FactoMineR::PCA(qData, 
                             scale.unit = var.scaling, 
                             ncp = ncp, 
                             graph= FALSE)
  
  return(res.pca)
}



#' @title Plots variables of PCA
#' 
#' @param res.pca Result of FactoMineR::PCA
#' 
#' @param chosen.axes The dimensions to plot
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
"plotPCA_Var"

#' @importFrom factoextra fviz_pca_var
#' 
#' @export
#' 
plotPCA_Var <- function(res.pca = NULL, 
                        chosen.axes = c(1,2)){
  if (is.null(res.pca))
    return(NULL)
  
  
  factoextra::fviz_pca_var(res.pca, 
                           axes = chosen.axes, 
                           col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE # Evite le chevauchement de texte
                           )
  
}



#' @title Plots individuals of PCA
#' 
#' @param res.pca Result of FactoMineR::PCA
#' 
#' @param chosen.axes The dimensions to plot
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
"plotPCA_Ind"

#' @importFrom factoextra fviz_pca_ind
#' 
#' @export
#' 
plotPCA_Ind <- function(res.pca, chosen.axes=c(1,2)){
  if (is.null(res.pca))
    return(NULL)
  
  factoextra::fviz_pca_ind(res.pca,  
                           axes = chosen.axes, 
                           geom  ="point")
  
}





#' @title Plots the eigen values of PCA with the highcharts library
#' 
#' @param res.pca Result of FactoMineR::PCA
#' 
#' @return A histogram
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
"plotPCA_Eigen_hc"

#' @export
#' 
#' @import highcharter
#' 
plotPCA_Eigen_hc <- function(res.pca){
  if (is.null(res.pca)){return(NULL)}
  hc <- highcharter::highchart() %>%
    hc_yAxis_multiples(list(title = list(text = "% of variances"),lineWidth = 0, labels = list(format = "{value}%"), max = 100), 
                       list(title = list(text = "Cumulative % of variances"), opposite = FALSE, max = 100),
                       list(title = list(text = "Eigen values"), opposite = TRUE, labels = list(format = "{value}%") )) %>%
    hc_xAxis(title = "Principal Components", categories = rownames(res.pca$eig)) %>%
    hc_add_series(data.frame(y=res.pca$eig[,2]), type="column",  name = "% of variances", yAxis = 0) %>%
    hc_add_series(data.frame(y=res.pca$eig[,3]), type="line",  color="darkblue",name = "Cumulative % of variances", marker = "diamond", color = "#FF7900", yAxis = 0) %>%
    #hc_add_series(data.frame(y=res.pca$eig[,1]),  type="line",  lineWidth = 4,name = "Eigen values", color = "orange", yAxis = 2) %>%
    #hc_tooltip(crosshairs = TRUE, headerFormat = "<b>{point.x}</b><br>") %>%
    hc_legend(enabled = TRUE)
  # hc_plotOptions(column = list(colorByPoint = TRUE, colors = SiteOTD$Colors))
  
}
