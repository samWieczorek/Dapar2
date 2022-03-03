#' 
#' #' This function computes few values about the adjacency matrix such as the number of proteins that are only defined by 
#' #' specific peptides, shared peptides or a mixture of two. 
#' #' 
#' #' @title Computes the number of proteins that are only defined by 
#' #' specific peptides, shared peptides or a mixture of two.
#' #' 
#' #' @param X The adjacency matrix with both specific and shared peptides.
#' #' 
#' #' @return A list of values:
#' #' * nbPeptides: the number of peptides in the matrix,
#' #' nbSpecificPeptides: the number of specific peptides in the matrix,
#' #' nbSharedPeptides: the number of shared peptides in the matrix,
#' #' nbProt: the number of proteins in the matrix,
#' #' protOnlyUniquePep: the list of proteins only defined by specific peptides,
#' #' protOnlySharedPep: the list of proteins only defined by shared peptides,
#' #' protMixPep: the list of proteins defined by both shared and specific peptides.
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' #' obj <- Exp1_R25_pept[seq_len(1000),]
#' #' obj <- addListAdjacencyMatrices(obj, 2)
#' #' X <- GetAdjMat(obj[[2]])$all
#' #' matAdjStats(X)
#' #' 
#' #' @export
#' #' 
#' matAdjStats <- function(X){
#'   if (is.null(X)){
#'     warning('The adjacency matrix is NULL.')
#'     return(NULL)
#'   }
#'   
#'   
#'   ind.shared.Pep <- which(rowSums(as.matrix(X))>1)
#'   ind.unique.Pep <- which(rowSums(as.matrix(X))==1)
#'   
#'   M.shared.Pep <- X[ind.shared.Pep,]
#'   M.shared.Pep <- M.shared.Pep[,-which(colSums(as.matrix(M.shared.Pep))==0)]
#'   
#'   M.unique.Pep <- X[ind.unique.Pep,]
#'   M.unique.Pep <- M.unique.Pep[,-which(colSums(as.matrix(M.unique.Pep))==0)]
#'   
#'   
#'   pep.names.shared <- colnames(M.shared.Pep)
#'   pep.names.unique <- colnames(M.unique.Pep)
#'   protOnlyShared <- setdiff(pep.names.shared, intersect(pep.names.shared, pep.names.unique))
#'   protOnlyUnique <- setdiff(pep.names.unique, intersect(pep.names.shared, pep.names.unique))
#'   protMix <- intersect(pep.names.shared, pep.names.unique)
#'   
#'   return (list(nbPeptides = nrow(M.unique.Pep)+nrow(M.shared.Pep),
#'                nbSpecificPeptides = nrow(M.unique.Pep),
#'                nbSharedPeptides = nrow(M.shared.Pep),
#'                nbProt = length(protOnlyShared)+length(protOnlyUnique)+length(protMix),
#'                protOnlyUniquePep =protOnlyUnique,
#'                protOnlySharedPep =protOnlyShared,
#'                protMixPep = protMix))
#' }





#' Method to plot the disrtibution (histogram) of peptides w.r.t the proteins with proteins and peptides 
#' in an adjacency matrix
#' 
#' @title Function to create a histogram that shows the repartition of
#' peptides w.r.t. the proteins
#' 
#' @param X An adjacency matrix.
#' 
#' @param type A string which is the type of matrix (used to build the plot title). Default value is 'all'.
#' 
#' @return A histogram  
#' 
#' @author Alexia Dorffer, Samuel Wieczorek
#' 
#' @examples
#' data(ft)
#' X <- adjacencyMatrix(ft[[1]])
#' GraphPepProt_hc(X)
#' 
#' @import highcharter
#' 
#' @export
#' 
GraphPepProt_hc <- function(X, type = 'all'){
  if (is.null(X)){
    warning("'X' is empty.")
    return (NULL)
  } 
  X <- as(X, 'matrix')
  t <- t(X)
  t <- apply(X, 2, sum, na.rm=TRUE)
  tab <- table(t)
  conds <- names(tab)
  
  h1 <-  highchart() %>%
    customChart(chartType = "column") %>%
    hc_title(text = paste0("Distribution of ",
                           type, 
                           " peptides w.r.t. proteins")) %>%
    hc_add_series(data = tab, 
                  type = "column", 
                  colorByPoint = TRUE) %>%
    hc_colors('orange') %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation = list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds, 
             title = list(
               text = "Number of peptides")
             ) %>%
    hc_yAxis(categories = conds, 
             title = list(text = "Number of proteins")
             ) %>%
    customExportMenu(fname = "HistoMatAdj") %>%
    hc_tooltip(headerFormat = '',
               pointFormat = "{point.y}")
  
  h1
  
}

