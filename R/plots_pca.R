

#' Compute the PCA
#' 
#' @title Compute the PCA
#' 
#' @param obj numeric matrix
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
#' library(FactoMineR)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' res.pca <- wrapper.pca(qData)
#' 
#' @importFrom FactoMineR PCA
#' @importFrom stats na.omit
#' 
#' @export
wrapper.pca <- function(qData, var.scaling=TRUE, ncp=NULL){
  
  if (is.null(var.scaling)) {var.scaling <- TRUE}
  if (length(which(is.na(qData))) > 0){ qData <- stats::na.omit(qData) }
  
  if (is.null(ncp)){
    nmax <- 12
    y <- qData
    nprot <- dim(y)[1]
    n <- dim(y)[2] # If too big, take the number of conditions.
    
    if (n > nmax){
      n <- length(unique(Biobase::pData(obj)$Condition))
    }
    
    ncp <- min(n, nmax)
  }
  
  res.pca <- FactoMineR::PCA(qData, scale.unit = var.scaling, ncp=ncp, graph=FALSE)
  # si warning pour les missing values, le reproduire dans l'interface graphique
  
  return(res.pca)
}



#' Plots the variables of PCA
#' 
#' @title Plots variables of PCA
#' @param res.pca Result of FactoMineR::PCA
#' @param chosen.axes The dimensions to plot
#' @return A plot
#' @author Samuel Wieczorek, Enora Fremy
#' @examples
#' library(FactoMineR)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' res.pca <- wrapper.pca(qData)
#' plotPCA_Var(res.pca)
#' @importFrom factoextra fviz_pca_var
#' @export
plotPCA_Var <- function(res.pca, chosen.axes=c(1,2)){
  #plot.PCA(res.pca, choix="var", axes = chosen.axes, title="Sample factor map (PCA)")
  #require(factoextra)
  # Colorer en fonction du cos2: qualit? de repr?sentation
  if (is.null(res.pca)){
    return(NULL)
  }
  
  factoextra::fviz_pca_var(res.pca, axes = chosen.axes, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE # ?vite le chevauchement de texte
  )
  
}


#' Plots the individuals of PCA
#' 
#' @title Plots individuals of PCA
#' @param res.pca Result of FactoMineR::PCA
#' @param chosen.axes The dimensions to plot
#' @return A plot
#' @author Samuel Wieczorek, Enora Fremy
#' @examples
#' library(FactoMineR)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' res.pca <- wrapper.pca(qData)
#' plotPCA_Ind(res.pca)
#' @importFrom factoextra fviz_pca_ind
#' @export
plotPCA_Ind <- function(res.pca, chosen.axes=c(1,2)){
  #plot.PCA(res.pca, choix="ind", axes = chosen.axes, select = 0:-1, title="Protein factor map (PCA)")
  if (is.null(res.pca)){return(NULL)}
  #require(factoextra)
  factoextra::fviz_pca_ind(res.pca,  axes = chosen.axes, geom="point")
  
}




#' Plots the eigen values of PCA with the highcharts library
#' 
#' @title Plots the eigen values of PCA with the highcharts library
#' @param res.pca Result of FactoMineR::PCA
#' @return A histogram
#' @author Samuel Wieczorek, Enora Fremy
#' @examples
#' library(FactoMineR)
#' library(highcharter)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' res.pca <- wrapper.pca(qData, ncp=6)
#' plotPCA_Eigen_hc(res.pca)
#' @export
#' @import highcharter
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
