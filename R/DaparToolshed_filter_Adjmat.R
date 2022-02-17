##' @title Filter a peptide assay on the basis of its adjacency matrix.
##' 
##' @description
##' 
##' These functions filters (delete) peptides of an assay, applying a function
##' on peptides and proteins. They can be used alone but the usual usage is to
##' create an instance of a class [FunctionFilter] and to pass it to the function
##' [filterFeaturesOneSE] in order to create a new assay, embedded into the [QFeatures]
##' object.
##' 
##' @param object An object of class `SummarizedExperiment`
##' 
##' @param fun A `list()` of additional parameters
##' 
##' @param top A `integer(1)` which is the number of xxx
##' 
##' @details
##' 
##' This function builds an intermediate matrix with scores for each peptide
##' based on 'fun' parameter. Once this matrix is built, one select the 'n' peptides
##' which have the higher score
##' 
##' The list of filter functions is given by [adjMatFilters()]:
##'  
##' - [specPeptides()]: returns a new assay of class [SummazizedExperiment] with only
##' specific peptides;
##' 
##' - [sharedpeptides()]: returns a new assay of class [SummazizedExperiment] with only
##' shared peptides;
##' 
##' - [topnPeptides()]: returns a new assay of class [SummazizedExperiment] with only
##' the 'n' peptides which best satisfies the condition. The condition is represented by
##' functions which calculates a score for each peptide among all samples. The list
##' of these functions is given by [topnFunctions()]:
##' 
##' - [rowMedians()]: xxx;
##' 
##' - [rowMeans()]: xxx;
##' 
##' - [rowSums()]: xxx;
##' 
##' 
##' @seealso The [QFeatures-filtering-oneSE] man page for the class `FunctionFilter`.
##' 
##' @rdname adjacency-matrix-filter
##'
##' @author Samuel Wieczorek
##' 
##' @examples 
##' 
##' #------------------------------------------------
##' # This function will keep only specific peptides
##' #------------------------------------------------
##' 
##' FunctionFilter('specPeptides', list()))
##' 
##' #------------------------------------------------
##' # This function will keep only shared peptides
##' #------------------------------------------------
##' 
##' FunctionFilter('sharedPeptides', list()))
##' 
##' #------------------------------------------------
##' # This function will keep only the 'n' best peptides
##' # w.r.t the quantitative sum of each peptides among 
##' # all samples
##' #------------------------------------------------
##' 
##' FunctionFilter('topnPeptides', fun = 'rowSums', top = 2)
##'
NULL


##' @export
##' @rdname adjacency-matrix-filter
adjMatFilters <- function()
  c('specPeptides', 'sharedPeptides', 'topnPeptides')



##' @export
##' @rdname adjacency-matrix-filter
specPeptides <- function(object){
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


##' @export
##' @rdname adjacency-matrix-filter
sharedPeptides <- function(object){
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


##' @export
##' @rdname adjacency-matrix-filter
topnFunctions <- function()
  c('rowMedians', 'rowMeans', 'rowSums')


##' @export
##' @rdname adjacency-matrix-filter
topnPeptides <- function(object, fun, top){
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(object)))
  
  # Preparing the variables
  #fun <- list(...)[[1]]$fun
  #n <- list(...)[[1]]$top
  X <- adjacencyMatrix(object)
  qData <- assay(object)
  
  stopifnot(inherits(X, "dgCMatrix"))
  
  if(!(fun %in% topnFunctions())){
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
