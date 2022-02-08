##' @exportMethod filterXFeatures
##' @rdname QFeatures-X-filtering
setMethod("filterXFeatures", "QFeatures",
          function(object, i, name = "newAssay", mode = 'all', ...) {
            if (isEmpty(object))
              return(object)
            if (name %in% names(object))
              stop("There's already an assay named '", name, "'.")
            if (missing(i))
              i <- main_assay(object)
            
            ## Create the aggregated assay
            filtered.SE <- filterXFeatures(object[[i]], mode, ...)
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               filtered.SE,
                               name = name)
            ## Link the input assay to the aggregated assay
            rowData(object[[i]])[['_temp_ID_']] <- rownames(object[[i]])
            rowData(object[[name]])[['_temp_ID_']] <- rownames(object[[name]])
            addAssayLink(object,
                         from = i,
                         to  = name,
                         varFrom = '_temp_ID_',
                         varTo = '_temp_ID_')
            
            rowData(object[[i]])[['_temp_ID_']] <- NULL
            rowData(object[[name]])[['_temp_ID_']] <- NULL
            
            
          })


##' @exportMethod aggregateFeatures
##' @rdname QFeatures-X-filtering
setMethod("filterXFeatures", "SummarizedExperiment",
          function(object, mode = 'all', ...){
            
            # Compute the theoretical new adjacency matrix
            X <- adjacencyMatrix(object)
            X.new <- updateAdjacencyMatrix(X, mode, ...)
            
            # Find and delete lonely peptides
            lonelyPeptides <- which(rowSums(X) == 0)
            object <- object[-lonelyPeptides]
            
              # Find and delete lonely proteins
            lonelyProteins <- which(colSums(X) == 0)
            names.lonelyProts <- colnames(X)[lonelyProteins]
            
            rd <- unfoldDataFrame(rowData(feat2[[2]]), fcol)
            ind <- match(names.lonelyProts, rd[, fcol])
            rd <- rd[-ind,]
            rd <- reduceDataFrame(rd, rd[,fcol], drop = TRUE)
            }
          )


.Get_Indices_Filter_X <- function(X, mode = 'all', ...){
  
  X.binary <- X
  # Used if X is a weighted matrix
  X.binary[which(X.binary != 0)] <- 1
  ind <- seq_len(nrow(X))
  ind <- switch(mode,
                spec = which(rowSums(X.binary) == 1),
                shared = which(rowSums(X.binary) > 1),
                topn = X <- .Build_Topn_Mat(X, ...)
  )
  return(ind)
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
.Build_Topn_Mat <- function(X, qData = NULL, fun = 'rowMedians', n = 10){
  
  stopifnot(inherits(X, "dgCMatrix"))
  if(!(fun %in% c('rowMedians', 'rowMeans', 'rowSums'))){
    warning("'fun' must be one of the following: 'rowMedians', 'rowMeans', 'rowSums'")
    return(NULL)
  }
  
  temp.X <- as(X * do.call(fun, list(qData)), "dgCMatrix")
  
  # Get the 'n' entities with the best score for each column
  for (c in seq_len(ncol(X))){
    v <- order(temp.X[,c],decreasing=TRUE)[seq_len(n)]
    l <- v[which((temp.X[,c])[v] != 0)]
    
    if (length(l) > 0){
      diff <- setdiff( which(X[,c] == 1), l)
      if (length(diff)) {X[diff,c] <- 0}
    }
  }
  
  X
}



