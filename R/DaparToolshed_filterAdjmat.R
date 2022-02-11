##' @exportMethod filterAdjmatOneSE
##' @rdname AdjMat-filtering
setMethod("filterAdjmatOneSE", "QFeatures",
          function(object, i, name = "newAssay", idcol, fcol, filters, ...) {
            if (isEmpty(object))
              return(object)
            if (name %in% names(object))
              stop("There's already an assay named '", name, "'.")
            if (missing(i))
              i <- main_assay(object)
            #browser()
            if (missing(idcol))
              idcol <- '_temp_ID_'
            else {
              if (!is.null(idcol)){
                if (nrow(object[[i]]) != length(unique(rowData(object[[i]])[, idcol]))){
                warning("'idcol' must contain different values. As it is not the case, 
                        the id of the dataset will be represented by the column
                        named '_temp_ID_'")
                idcol <- '_temp_ID_'
                }
              }
            }
            ## Create the aggregated assay
            filtered.SE <- filterAdjmatOneSE(object[[i]], fcol, filters, ...)
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               filtered.SE,
                               name = name)
            ## Link the input assay to the aggregated assay
            rowData(object[[i]])[idcol] <- rownames(object[[i]])
            rowData(object[[name]])[idcol] <- rownames(object[[name]])
            object <- addAssayLink(object,
                         from = names(object)[i],
                         to  = name,
                         varFrom = idcol,
                         varTo = idcol)
            
            #browser()
            #rowData(object[[i]])[['_temp_ID_']] <- NULL
            #rowData(object[[name]])[['_temp_ID_']] <- NULL
            
           return(object) 
          })


##' @exportMethod filterAdjmatOneSE
##' @rdname AdjMat-filtering
setMethod("filterAdjmatOneSE", "SummarizedExperiment",
          function(object, fcol, filters, ...){
             stopifnot('adjacencyMatrix' %in% names(rowData(object)))
              # Reduce the adjacency matrix w.r.t the modes parameters
              # Apply a filter at a time and update the object
              for (p in names(filters)){
                X <- adjacencyMatrix(object)
                filtered.X <- updateAdjacencyMatrix(X = X, mode = p, filters[[p]])
                rowData(object)$adjacencyMatrix <- filtered.X
                
                lonelyPeptides <- which(rowSums(as.matrix(filtered.X)) == 0)
                if (length(lonelyPeptides) > 0){
                  object <- object[-lonelyPeptides]
                  filtered.X <- adjacencyMatrix(object)
                }
                
                # Find and delete lonely proteins
                
                lonelyProteins <- which(colSums(as.matrix(filtered.X)) == 0)
                if (length(lonelyProteins) > 0){
                  X <- filtered.X[, -lonelyProteins]
                  rowData(object)$adjacencyMatrix <- X
                  rowData(object)[,fcol] <- makePeptideProteinVector(X)
                }
              }
              
              
              
              return(object)

            }
          )




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
Build_Topn_Mat <- function(X, qData = NULL, fun = 'rowMedians', n = 10){
  
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

  return(X)
}




#' @title xxx
#' 
#' @description xxx
#' @details Mode can be
#' - all xxxx
#' - onlySpec xxx
#' - onlyShared xxx
#' - topn xxx
#' 
#' @param se xxx
#' @param X xxx
#' @param mode xxx
#' @param ... Additional parameters passed to some mode functions
#' 
#' @export
#' 
#' @examples
#' feat2 <- readRDS('~/GitHub/DaparToolshedData/data/Exp2_R100_pept.rda')
#' feat2 <- feat2[1:10,]
#' X <- makeAdjacencyMatrix(rowData(feat2[[2]])[,'Protein_group_IDs'])
#' rownames(X) <- rownames(feat2[[2]])
#' updateAdjacencyMatrix(X, mode = 'all')
#'
updateAdjacencyMatrix <- function(X, mode, ...){
  
  argg <- list(...)[[1]]
  X.binary <- as.matrix(X) != 0
  # Used if X is a weighted matrix
  switch(mode,
         spec = {
           ind <- which(rowSums(X.binary) > 1)
           if (length(ind) > 0)
             X[ind,] <- 0
           },
         shared = {
           ind <- which(rowSums(X.binary) == 1)
           if (length(ind) > 0)
             X[ind,] <- 0
           },
         topn = {
           X <- Build_Topn_Mat(X = X,
                               qData = argg[['qData']],
                               fun = argg[['fun']],
                               n = argg[['n']])
         }
  )
  return(X)
}


# 
# 
# ##' @exportMethod aggregateFeatures
# ##' @rdname QFeatures-filtering-oneSE
# setMethod("filterFeaturesOneSE", "QFeatures",
#           function(object, i, name, filters.def, ...) {
#             if (isEmpty(object))
#               return(object)
#             if (name %in% names(object))
#               stop("There's already an assay named '", name, "'.")
#             if (missing(i))
#               i <- main_assay(object)
#             
#             ## Create the aggregated assay
#             aggAssay <- filterFeaturesOneSE(object[[i]], fcol, fun, ...)
#             
#             new.object <- filterFeatures()
#             ## Add the assay to the QFeatures object
#             object <- addAssay(object,
#                                new.object[length(new.object)],
#                                name = name)
#             ## Link the input assay to the aggregated assay
#             addAssayLink(object,
#                          from = i,
#                          to  = name,
#                          varFrom = fcol,
#                          varTo = fcol)
#           })
# 

