
##' @exportClass NumericVariableFilter
##' @rdname QFeatures-filtering
setClass("ComplexFilter",
         slots = c(name="character", params="list"),
         prototype = list(name=character(), params=list())
)


##' @param field `character(1)` refering to the name of the variable
##'     to apply the filter on.
##'
##' @param value `character()` or `integer()` value for the
##'     `CharacterVariableFilter` and `NumericVariableFilter` filters
##'     respectively.
##'
##' @param condition `character(1)` defining the condition to be used in
##'     the filter. For `NumericVariableFilter`, one of `"=="`,
##'     `"!="`, `">"`, `"<"`, `">="` or `"<="`. For
##'     `CharacterVariableFilter`, one of `"=="`, `"!="`,
##'     `"startsWith"`, `"endsWith"` or `"contains"`. Default
##'     condition is `"=="`.
##'
##' @param not `logical(1)` indicating whether the filtering should be negated
##'     or not. `TRUE` indicates is negated (!). `FALSE` indicates not negated.
##'     Default `not` is `FALSE`, so no negation.
##'
##' @export VariableFilter
##' @rdname QFeatures-filtering
ComplexFilter <- function(name, params) {
  new("ComplexFilter",
        name = name,
        params = params)
}




##' @exportMethod filterFeaturesOneSE
##' @rdname AdjMat-filtering
setMethod("filterFeaturesOneSE", "QFeatures",
          function(object, i, name = "newAssay", filters) {
            if (isEmpty(object))
              return(object)
            if (name %in% names(object))
              stop("There's already an assay named '", name, "'.")
            if (missing(i))
              i <- main_assay(object)
            
            if (missing(filters))
              return(object)
            
            if (is.null(metadata(object[[i]])$idcol)){
              warning('xxx')
              metadata(object[[i]])$idcol <- '_temp_ID_'
            }

            ## Create the aggregated assay
            new.se <- filterFeaturesOneSE(object[[i]], filters)
           
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               new.se,
                               name = name)
            
            if (nrow(new.se) > 0){
              idcol <- metadata(object[[i]])$idcol
              ## Link the input assay to the aggregated assay
              rowData(object[[i]])[,idcol] <- rownames(object[[i]])
              rowData(object[[name]])[,idcol] <- rownames(object[[name]])
              object <- addAssayLink(object,
                         from = names(object)[i],
                         to  = name,
                         varFrom = idcol,
                         varTo = idcol)
            }
            
           return(object) 
          })


##' @exportMethod filterFeaturesOneSE
##' @rdname AdjMat-filtering
setMethod("filterFeaturesOneSE", "SummarizedExperiment",
          function(object, filters){
            for (f in filters)
              object <- do.call(f@name, list(object, f@params))
            return(object)
            }
          )


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
  emptyLines <- which(rowSums(X) == 0)
  if (length(emptyLines) > 0)
    object <- object[-emptyLines]

  # Reload the adjacency matrix after lines deletion
  # Identify empty columns in the adjacency matrix
  X <- adjacencyMatrix(object)
  emptyCols <- which(colSums(X) == 0)

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
