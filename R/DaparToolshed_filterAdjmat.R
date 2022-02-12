

filtersOnAdjmat <- function()
  c('specPeptides', 'sharedPeptides', 'topnPeptides')

specPeptides <- function(object, ...){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(object)))
  stopifnot(!is.null(metadata(object)$fcol))
  
  X <- adjacencyMatrix(object)
  
  # Mask for only specific peptides
  x.spec <- as.matrix(X) != 0
  ind <- which(rowSums(x.spec) > 1)
  if (length(ind) > 0){
    X[ind,] <- 0
    object <- .UpdateSEBasedOnAdjmat(object, X)
  }
  
  return(object)
}


sharedPeptides <- function(object, ...){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(object)))
  stopifnot(!is.null(metadata(object)$fcol))
  
  X <- adjacencyMatrix(object)
  
  # Mask for only shared peptides
  x.spec <- as.matrix(X) != 0
  ind <- which(rowSums(x.spec) == 1)
  if (length(ind) > 0){
    X[ind,] <- 0
    object <- .UpdateSEBasedOnAdjmat(object, X)
  }
  
  return(object)
}



.UpdateSEBasedOnAdjmat <- function(object, X){
  rowData(object)$adjacencyMatrix <- X
  # Identify and delete the empty lines in the dataset
  emptyLines <- which(rowSums(as.matrix(X)) == 0)
  if (length(emptyLines) > 0)
    object <- object[-emptyLines]

  # Reload the adjacency matrix after lines deletion
  # Identify empty columns in the adjacency matrix
  X <- adjacencyMatrix(object)
  emptyCols <- which(colSums(as.matrix(X)) == 0)

  if (length(emptyCols) > 0){
    X <- X[, -emptyCols]
    rowData(object)$adjacencyMatrix <- X
    rowData(object)[,metadata(object)$fcol] <- makePeptideProteinVector(X)
  }

return(object)
}


#' @title xxxxx
#' @description xxx 
#' @details This function builds an intermediate matrix with scores for each peptide
#' based on 'fun' parameter. Once this matrix is built, one select the 'n' peptides
#' which have the higher score
#' 
#' - rowMedians xxx
#' - rowMeans xxx
#' - rowSums xxx
#' 
#' @param qData xxx
#' @param X xxx
#' @param fun xxx
#' @param n xxx
#' 
#' 
#' @examples 
#'
#' @export
#' 
topnPeptides <- function(object, ...){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(object)))
  
  # Preparing the variables
  fun <- list(...)[[1]]$fun
  n <- list(...)[[1]]$n
  X <- adjacencyMatrix(object)
  qData <- assay(object)
  
  stopifnot(inherits(X, "dgCMatrix"))
  
  if(!(fun %in% c('rowMedians', 'rowMeans', 'rowSums'))){
    warning("'fun' must be one of the following: 'rowMedians', 'rowMeans', 'rowSums'")
    return(NULL)
  }
  
  temp.X <- as(X * do.call(fun, list(qData)), "dgCMatrix")
  
  # Get the 'n' entities with the best score for each column
  for (c in seq_len(ncol(temp.X))){
    v <- order(temp.X[,c],decreasing=TRUE)[seq_len(n)]
    l <- v[which((temp.X[,c])[v] != 0)]
    if (length(l) > 0)
      X[-l, c] <- 0
  }
  object <- .UpdateSEBasedOnAdjmat(object, X)
  return(object)
}
