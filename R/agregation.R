#' 
#' 
#' 
#' 
#' #' Method to create a binary matrix with proteins in columns and peptides 
#' #' in lines on a \code{MSnSet} object (peptides)
#' #' 
#' #' @title Function matrix of appartenance group
#' #' @param pg A vector of xxxxx
#' #' @param names The names of peptides xxxx
#' #' @param unique A boolean to indicate whether only the unique peptides must 
#' #' be considered (TRUE) or if the shared peptides have to 
#' #' be integrated (FALSE).
#' #' @return A binary matrix  
#' #' @author Samuel Wieczorek
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' #' PG <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$parentProtId]
#' #' names <- names((Exp1_R25_pept[['original']]))
#' #' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' #' @export
#' BuildAdjacencyMatrix <- function(pg, names, unique=FALSE){
#'   PG.l <- strsplit(as.character(pg), split=";", fixed=TRUE)
#'   
#'   t <- table(data.frame(A=rep(seq_along(PG.l), lengths(PG.l)), B=unlist(PG.l)))
#'   
#'   if (unique == TRUE){
#'     ll <- which(rowSums(t)>1)
#'     if (length(ll) > 0) {
#'       t[ll,] <- 0
#'     }
#'   }
#'   
#'   X <- Matrix::Matrix(t, dimnames = list(names, colnames(t))
#'   )
#'   
#'   return(X)
#' }



#' Method to create a binary matrix with proteins in columns and peptides 
#' in lines on a \code{MSnSet} object (peptides)
#' 
#' @title Function matrix of appartenance group
#' @param pg A vector of xxxxx
#' @param names The names of peptides xxxx
#' @param type  A vector of type of adjacency matrices to include in the return value. Values are : 'All' for 
#' the matrix containing all peptides, 'Shared' for the matrix which contains only the shared peptides v=between proteins,
#' 'Specific' for the matrix which contains only the specific peptides. Default value for type is all three matrices
#' @return A list of three adjacency matrices. 1 - the matrix contains all peptides (shared and specific), 2 - 
#' the matrix contains only the peptides shared between proteins, 3 - the matrix contains only specific peptides 
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original']]))
#' X <- BuildListAdjacencyMatrices(PG, names, type='All')
#' @export
BuildListAdjacencyMatrices <- function(pg, names, type = c('All', 'Shared', 'Specific')){
  Xshared <- Xspec <- Xall <- NULL
  
  # the separator that is allowed is ','
  pg <- gsub(";", ",", as.character(pg), fixed=TRUE)
  PG.l <- strsplit(as.character(pg), split=",", fixed=TRUE)
  t <- table(data.frame(A=rep(seq_along(PG.l), lengths(PG.l)), B=unlist(PG.l)))
  
  #compute the matrix with all peptides
  if (sum(type=='All') == 1){
    Xall <- Matrix::Matrix(t, sparse=T,dimnames = list(names, colnames(t)))
  }
  
  #compute the matrix with only shared peptides
  if (sum(type=='Shared') == 1){
    tmpShared <- t
    ll <- which(rowSums(tmpShared)==1)
    if (length(ll) > 0) {
      tmpShared[ll,] <- 0
    }
    Xshared <- Matrix::Matrix(tmpShared, sparse=T, dimnames = list(names, colnames(t)))
  }
  
  #compute the matrix with only specific peptides
  if (sum(type=='Specific') == 1){
    tmpSpec <- t
    ll <- which(rowSums(tmpSpec)>1)
    if (length(ll) > 0) {
      tmpSpec[ll,] <- 0
    }
    Xspec <- Matrix::Matrix(tmpSpec, sparse=T, dimnames = list(names, colnames(t)))
  }
  
  return (list(all = Xall, onlyShared = Xshared, onlySpec = Xspec))
}






#' This function computes the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two. 
#' 
#' @title Computes the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two.
#' @param matShared The adjacency matrix with both specific and 
#' shared peptides.
#' @return A list
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' matAdjStats(X)
#' @export
matAdjStats <- function(matShared){
  if (is.null(matShared)){
    warning('The adjacency matrix is NULL.')
    return(NULL)
  }
  
  
  ind.shared.Pep <- which(rowSums(as.matrix(matShared))>1)
  ind.unique.Pep <- which(rowSums(as.matrix(matShared))==1)
  
  M.shared.Pep <- matShared[ind.shared.Pep,]
  M.shared.Pep <- M.shared.Pep[,-which(colSums(as.matrix(M.shared.Pep))==0)]
  
  M.unique.Pep <- matShared[ind.unique.Pep,]
  M.unique.Pep <- M.unique.Pep[,-which(colSums(as.matrix(M.unique.Pep))==0)]
  
  
  pep.names.shared <- colnames(M.shared.Pep)
  pep.names.unique <- colnames(M.unique.Pep)
  protOnlyShared <- setdiff(pep.names.shared, intersect(pep.names.shared, pep.names.unique))
  protOnlyUnique <- setdiff(pep.names.unique, intersect(pep.names.shared, pep.names.unique))
  protMix <- intersect(pep.names.shared, pep.names.unique)
  
  return (list(nbPeptides = nrow(M.unique.Pep)+nrow(M.shared.Pep),
               nbSpecificPeptides = nrow(M.unique.Pep),
               nbSharedPeptides = nrow(M.shared.Pep),
               nbProt = length(protOnlyShared)+length(protOnlyUnique)+length(protMix),
               protOnlyUniquePep =protOnlyUnique,
               protOnlySharedPep =protOnlyShared,
               protMixPep = protMix))
}




#' This function computes the number of peptides used to aggregate proteins.
#' 
#' @title Compute the number of peptides used to aggregate proteins
#' @param M A "valued" adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' @return A vector of boolean which is the adjacency matrix 
#' but with NA values if they exist in the intensity matrix.
#' @author Alexia Dorffer
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' CountPep(X)
#' @export
CountPep <- function (M) {
  z <- M
  z[z!=0] <- 1
  return(z)
}



#' Method to plot the disrtibution of peptides w.r.t the proteins with proteins and peptides on
#' a MSnSet object (peptides)
#' 
#' @title Function to create a histogram that shows the repartition of
#' peptides w.r.t. the proteins
#' @param mat An adjacency matrix.
#' @return A histogram  
#' @author Alexia Dorffer, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' GraphPepProt(X)
#' @export
GraphPepProt <- function(mat){
  if (is.null(mat)){
    warning('The parameter mat is empty.')
    return (NULL)
  } 
  
  #mat <- as.matrix(mat)
  t <- t(mat)
  t <- apply(mat, 2, sum, na.rm=TRUE)
  tab <- table(t)
  position <- seq(1, length(tab),by=3)
  conds <- names(tab)
  
  #par(mar=c(6,4,4,8) + 0.1)#, mgp=c(3,0.5,0)
  barplot(tab, 
          xlim=c(1, length(tab)),
          xlab="Nb of peptides", 
          ylab="Nb of proteins",
          names.arg=conds, 
          xaxp=c(1, length(tab), 3), 
          las=1
          , col = "orange")
  
}



#' Method to compute the detailed number of quantified peptides used for aggregating each protein
#' 
#' @title Computes the detailed number of peptides used for aggregating each protein 
#' @param X.list An adjacency matrix
#' @param qPepData A data.frame of quantitative data
#' @return A list of two items
#' @author Samuel Wieczorek
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qPepData <- assay(Exp1_R25_pept,'original_log')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X.list <- BuildListAdjacencyMatrices(PG, names, type=c('All', 'Shared', 'Specific'))
#' n <- GetDetailedNbPeptidesUsed(X.list, qPepData)
#' @export
GetDetailedNbPeptidesUsed <- function(X.list, qPepData){
  res <- list()
  
  qPepData[!is.na(qPepData)] <- 1
  qPepData[is.na(qPepData)] <- 0
  
  
  if (!is.null(X.list$all))
    res$nAll <- t(X.list$all) %*% qPepData
  if (!is.null(X.list$all))
    res$nShared <- t(X.list$onlyShared) %*% qPepData
  if (!is.null(X.list$all))
    res$nSpec <- t(X.list$onlySpec) %*% qPepData
  
  return(res)
}


#' Method to compute the detailed number of quantified peptides for each protein
#' 
#' @title Computes the detailed number of peptides for each protein 
#' @param X.list An adjacency matrix
#' @return A data.frame
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildListAdjacencyMatrices(PG, names, type=c('All','Shared', 'Specific'))
#' n <- GetDetailedNbPeptides(X)
#' @export
GetDetailedNbPeptides <- function(X.list){
  return(list(nTotal = rowSums(t(as.matrix(X.list$all))),
              nShared=rowSums(t(X.list$onlyShared)), 
              nSpec=rowSums(t(X.list$onlySpec)))
  )
  
}





#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData A data.frame of quantitative data
#' @param X An adjacency matrix
#' @return A matrix
#' @author Samuel Wieczorek
#' @example
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildListAdjacencyMatrices(PG, names, type=c('All'))
#' n <- inner.sum(assay(Exp1_R25_pept[['original_log']], X$all))
#' @export
inner.sum <- function(qPepData, X, na.rm=TRUE){
  qPepData[is.na(qPepData)] <- 0
  Mp <- t(X) %*% qPepData
  return(Mp)
}



#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData A data.frame of quantitative data
#' @param X An adjacency matrix
#' @return xxxxx
#' @author Samuel Wieczorek
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildListAdjacencyMatrices(PG, names, type=c('All'))
#' inner.mean(assay(Exp1_R25_pept[['original_log']], X$all)
#' @export
inner.mean <- function(qPepData, X){
  Mp <- inner.sum(qPepData, X)
  Mp <- Mp / GetNbPeptidesUsed(X, qPepData)
  
  return(Mp)
}





#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData A data.frame of quantitative data
#' @param X An adjacency matrix
#' @param method xxxxx
#' @param n xxxxx
#' @return xxxxx
#' @author Samuel Wieczorek
#' @export
#' @example 
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' qPepData <- assay(Exp1_R25_pept[['original_log']])
#' X <- BuildListAdjacencyMatrices(PG, names, type=c('All'))
#' inner.aggregate.topn(qPepData, X$all)
#' @importFrom stats median
inner.aggregate.topn <-function(qPepData, X, method='Mean', n=10){
  
  med <- apply(qPepData, 1, median)
  xmed <- as(X * med, "dgCMatrix")
  for (c in 1:ncol(X)){
    v <- order(xmed[,c],decreasing=TRUE)[1:n]
    l <- v[which((xmed[,c])[v] != 0)]
    
    if (length(l) > 0){
      diff <- setdiff( which(X[,c] == 1), l)
      if (length(diff)) {X[diff,c] <- 0}
    }
  }
  
  Mp <- NULL
  switch(method,
         Mean = Mp <- inner.mean(qPepData, X),
         Sum = Mp <- inner.sum(qPepData, X)
  )
  
  return(Mp)
}






#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData xxxxx
#' @param X xxxx
#' @param init.method xxx
#' @param method xxx
#' @param n xxxx
#' @return xxxxx
#' @author Samuel Wieczorek
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' qPepData <- assay(Exp1_R25_pept[['original_log']])
#' X <- BuildListAdjacencyMatrices(PG, names, type=c('All'))
#' inner.aggregate.iter(qPepData, X$all)
#' @export
inner.aggregate.iter <- function(qPepData, X, init.method='Sum', method='Mean', n=NULL){
  
  if (!(init.method %in% c("Sum", "Mean"))) {
    warning("Wrong parameter init.method")
    return(NULL)
  }
  
  if (!(method %in% c("onlyN", "Mean"))){
    warning("Wrong parameter method")
    return(NULL)
  }
  
  
  if (method=='onlyN' && is.null(n)){
    warning("Parameter n is null")
    return(NULL)
  }
  
  yprot <- NULL
  switch(init.method,
         Sum= yprot <- inner.sum(qPepData, X),
         Mean= yprot <- inner.mean(qPepData, X)
  )
  conv <- 1
  
  while(conv > 10**(-10)){
    mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot[is.na(mean.prot)] <- 0
    
    X.tmp <- mean.prot*X
    X.new <- X.tmp/rowSums(as.matrix(X.tmp), na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    # l'appel ? la fonction ci-dessous d?pend des param?tres choisis par l'utilisateur
    switch(method,
           Mean = yprot <- inner.mean(qPepData, X.new),
           onlyN = yprot <- inner.aggregate.topn(qPepData,X.new,'Mean', n)
    )
    
    mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot.new[is.na(mean.prot.new)] <- 0
    
    conv <- mean(abs(mean.prot.new - mean.prot))
    print(paste0("conv : ", conv))
  }
  return(as.matrix(yprot))
}




#' This function computes the intensity of proteins based on the sum of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins with the sum of the intensities
#' of their peptides.
#' @param qPepData A matrix of intensities of peptides
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' @return A matrix of intensities of proteins
#' @author Alexia Dorffer, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' aggSum(assay(Exp1_R25_pept,2), X)
#' @export
aggSum <- function(qPepData, X, ...){
  #qPepData <- 2^(qPepData)
  protData <- inner.sum(qPepData, X,...)
  #obj.prot <- finalizeAggregation(obj.pep, pepData, protData, X)
  return(protData)
}



#' This function computes the intensity of proteins as the mean of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins as the mean of the intensities
#' of their peptides.
#' @param qPepData A matrix of intensities of peptides
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' @return A matrix of intensities of proteins
#' @author Alexia Dorffer
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' aggMean(assay(Exp1_R25_pept,2), X)
#' @export
aggMean <- function(qPepData, X){
  #qPpepData <- 2^(qPepData)
  protData <- inner.mean(qPpepData, X)
  #obj.prot <- finalizeAggregation(obj.pep, pepData, protData, X)
  return(protData)
}




#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData xxxxx
#' @param X xxxx
#' @param conditions xxxx
#' @param init.method xxxxx
#' @param method xxxxx
#' @param n xxxx
#' @return xxxxx
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']][1:1000,])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']][1:1000]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' conditions <- colData(Exp1_R25_pept)@listData$Condition
#' aggIterParallel(assay(Exp1_R25_pept,2)[1:1000,], X, conditions)
#' @export
#' @importFrom doParallel registerDoParallel
#' @import foreach
aggIterParallel <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
  if (is.null(conditions)){
    warning('The parameter conds is NULL: the aggregation cannot be process.')
    return(NULL)
  }
  require(doParallel)
  doParallel::registerDoParallel()
  
  #qPepData <- 2^(qPepData)
  protData <- matrix(rep(0,ncol(X)*nrow(X)), nrow=ncol(X))
  
  protData <- foreach(cond = 1:length(unique(conditions)), .combine=cbind) %dopar% {
    condsIndices <- which(conditions == unique(conditions)[cond])
    qData <- qPepData[,condsIndices]
    inner.aggregate.iter(qData, X, init.method, method, n)
  }
  
  protData <- protData[,colnames(qPepData)]
  #obj.prot <- finalizeAggregation(obj.pep, qData.pep, protData, X)
  
  return(protData)
  
}




#' Method to xxxxx
#' 
#' @title xxxx 
#' @param qPepData xxxxx
#' @param X xxxx
#' @param conditions xxx
#' @param init.method xxxxx
#' @param method xxxxx
#' @param n xxxx
#' @return A protein object of class \code{MSnset}
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']][1:1000,])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']][1:1000]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' conditions <- colData(Exp1_R25_pept)@listData$Condition
#' aggIter(assay(Exp1_R25_pept,2)[1:1000,], X, conditions)
#' @export
aggIter <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
  if (is.null(conditions)){
    warning('The parameter conds is NULL: the aggregation cannot be process.')
    return(NULL)
  }
  
  ### a reproduire iterativement pour chaque condition
  # Initialisation: presque aucune d?pendance ? l'initialisation prendre "sum overall" et  matAdj = X par simplicit?
  #X <- as.matrix(X)
  #qPepData <- 2^(qPepData)
  
  protData <- matrix(rep(0,ncol(X)*ncol(qPepData)), nrow=ncol(X))
  
  for (cond in unique(conditions)){
    condsIndices <- which(conditions == cond)
    qData <- qPepData[,condsIndices]
    protData[,condsIndices]  <- inner.aggregate.iter(qData, X, init.method, method, n)
  }
  #obj.prot <- finalizeAggregation(obj.pep, qData.pep, protData, X)
  return(protData)

}




#' This function computes the intensity of proteins as the sum of the 
#' intensities of their n best peptides.
#' 
#' @title Compute the intensity of proteins as the sum of the 
#' intensities of their n best peptides.
#' @param qPepData A matrix of intensities of peptides
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' @param method xxx
#' @param n The maximum number of peptides used to aggregate a protein.
#' @return A matrix of intensities of proteins
#' @author Alexia Dorffer, Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']][1:1000,])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']][1:1000]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' aggTopn(assay(Exp1_R25_pept,2)[1:1000,], X, n=3)
#' @export
aggTopn <- function(qPepData, X,  method='Mean', n=10){
  #qPepData <- 2^(qPepData)
  
  protData <- inner.aggregate.topn(qPepData, X, method=method, n)
  
  #obj.prot <- finalizeAggregation(obj.pep, pepData, protData, X)
  return(protData)
}





##' @title Aggreagate quantitative features.
##' 
##' @description
##' This function takes a matrix of quantitative features `x` and a
##' factor (of length equal to `nrow(x)`) defining subsets, and
##' applies a user-defined function to aggregate each subset into a
##' vector of quantitative values. This function is the same as aggregate_by_vector
##'
##' User-defined functions must thus return a vector of length equal
##' to `ncol(x)`. Examples thereof are
##'
##' - [medianPolish()] to fits an additive model (two way decomposition)
##'   using Tukey's median polish_ procedure using
##'   [stats::medpolish()];
##'
##' - [robustSummary()] to calculate a robust aggregation using
##'   [MASS::rlm()];
##'
##' - [base::colMeans()] to use the mean of each column;
##'
##' - [base::colSums()] to use the sum of each column;
##'
##' - [matrixStats::colMedians()] to use the median of each column.
##'
##' @param qPepData xxxx. 
##' @param X Axxxxx.
##' @param FUN A `function` to be applied to the subsets of `x`.
##' @param ... Additional arguments passed to `FUN`.
##' @return A new `matrix` of dimensions `ncol(x)` and `length(INDEX)`
##'     with `dimnames` equal to `colnames(x)` and `INDEX`.
##' 
##' @author Samuel Wieczorek
##'
##' @family Quantitative feature aggregation
##' 
##' @export
aggregate_with_matAdj <- function(qPepData, X, FUN, ...){
  res <- do.call(FUN, list(qPepData, X, ...))
  res
}




##' @title Aggregate an assay's quantitative features which take into account
##' the peptides shared between proteins
##'
##' @description
##' 
##' This function aggregates the quantitative features of an assay,
##' applying an aggregation function (`fun`) to sets of features as
##' defined by the `fcol` feature variable. The new assay's features
##' will be named based on the unique `fcol` values.
##' This function is largely inspired by xxxx . The difference is that it provides
##' a mean to take into account the peptides shared between proteins.
##'
##'
##' @param object An instance of class [Features].
##'
##' @param i The index or name of the assay which features will be
##'     aggregated the create the new assay.
##'
##' @param X.list An adjacency matrix as computed by xxxxx.
##'
##' @param typeMatAdj xxx
##' 
##' @param name A `character(1)` naming the new assay. Default is
##'     `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##'     
##' @param meta.names A vector of character strings that are the metadata of the peptides which needs to be aggregated
##' and kept in the protein dataset
##'
##' @param fun A function used for quantitative feature
##'     aggregation. See Details for examples.
##'
##' @param ... Additional parameters passed the `fun`.
##'
##' @return A `Features` object with an additional assay.
##'
##' @details
##'
##' Aggregation is performed by a function that takes a matrix as
##' input and returns a xxxxx. Examples
##' thereof are
##'
##' - [DAPAR2:aggSum()] to use the sum of each column (default);
##'
##' - [DAPAR2:aggMean()] to use the sum of each column;
##'
##' - [DAPAR2:aggIter()] to use the mean of each column;
##'
##' - [DAPAR2:aggIterParallel()] same as previous function but use parallelism.
##'
##' - [DAPAR::aggTopn] to use the sum of each column;
##'
#' 
#' @seealso The *Features* vignette provides an extended example and
#'     the *Processing* vignette, for a complete quantitative
#'     proteomics data processing pipeline.
#' 
#' @aliases aggregateFeatures aggregateFeatures,Features-method
#'
#' @name aggregateFeatures
#'
#' @rdname Features-aggregate_dapar
#'
#' @importFrom MsCoreUtils aggregate_by_vector robustSummary
#'
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X <- BuildAdjacencyMatrix(PG, names, TRUE)
#' X.list <- BuildListAdjacencyMatrices(PG, names, type=c('All','Shared', 'Specific'))
#' conditions <- colData(Exp1_R25_pept)@listData$Condition
#' aggregateFeatures_sam(Exp1_R25_pept,2, X.list, typeMatAdj = 'All', name='aggregated', meta.names = 'Sequence', aggSum)
#' @export aggregateFeatures_sam
##-------------------------------------------------------------------------
# setMethod("aggregateFeatures_sam", "Features",
#           function(object, i, X, name = "newAssay",
#                    fun = aggSum, ...)
#             .aggregateFeatures_sam(object, i, fcol, name, fun, ...))


aggregateFeatures_sam <- function(object, i, X.list, typeMatAdj, name, meta.names = NULL, fun, ...) {
  if (isEmpty(object))
    return(object)
  if (name %in% names(object))
    stop("There's already an assay named '", name, "'.")
  if (missing(X.list))
    stop("'X.list' is required.")    
  if (missing(i))
    i <- main_assay(object)
  
  assay_i <- assay(object, i)
  rowdata_i <- rowData(object[[i]])

  
  ## Message about NA values is quant/row data
  has_na <- character()
  if (anyNA(assay_i))
    has_na <- c(has_na, "quantitative")
  if (anyNA(rowdata_i, recursive = TRUE))
    has_na <- c(has_na, "row")
  if (length(has_na)) {
    msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                       " data contain missing values."),
                 "Please read the relevant section(s) in the",
                 "aggregateFeatures manual page regarding the",
                 "effects of missing values on data aggregation.")
    message(paste(strwrap(msg), collapse = "\n"))
  }
  
  #aggregated_assay <- aggregate_by_vector(assay_i, groupBy, fun, ...)
  X <- NULL
  switch(typeMatAdj,
        All = X <- X.list$all,
        Shared = X <- X.list$onlyShared,
        Specific = X <- X.list$onlySpecific
        )
  aggregated_assay <- aggregate_with_matAdj(assay_i, X, fun, ...)
  
  
  # aggregated_rowdata <- Features::reduceDataFrame(rowdata_i, rowdata_i[[fcol]],
  #                                                 simplify = TRUE, drop = TRUE,
  #                                                 count = TRUE)
  ## correspond aux fData des proteines nouvellement creees
  aggregated_rowdata <- rowdata_stats_Aggregation_sam(assay_i, X.list)
  if (!is.null(meta.names)){
  aggregated_rowdata <- cbind(aggregated_rowdata, aggMetadata_sam(rowdata_i, meta.names, X))
  }

  # 
   se <- SummarizedExperiment(aggregated_assay,
                              rowData = aggregated_rowdata[rownames(aggregated_assay), ])
  
  # hits <- findMatches(rownames(aggregated_assay), groupBy)
   # rownames(aggregated_assay) : correspond au nom des proteines nouvellement creees
   # groupBy : correspond aaux id des peptides d'origine ?
   #hits <- findMatches(rownames(aggregated_assay), groupBy)
   ## from et to forment la definition du graphe pepetide-protein
   test <- which(as.matrix(X.list$all)==1, arr.ind=TRUE)
   
   from <- test[,'col']
   to <- test[,'row']
   hits <- Hits(from=from, to=to, nLnode=length(from), nRnode=nrow(X.list$all),sort.by.query=TRUE)
   
   elementMetadata(hits)$names_from <- rownames(assay_i)[hits@to]
   elementMetadata(hits)$names_to <- colnames(X.list$all)[hits@from]
  
   
   assayLinks <- AssayLink(name = name,
                           from = ifelse(is.character(i), i, names(object)[i]),
                           fcol = fcol,
                           hits = hits)
   addAssay(object,
            se,
            name = name,
            assayLinks = assayLinks)
   
}






#' Method to finalize the aggregation process
#' 
#' @title Finalizes the aggregation process 
#' @param qPepData A data.frame of quantitative data not logged of peptides
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' @return A protein object of class \code{SummarizedExperiment}
#' @author Samuel Wieczorek
#' @export
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X.list <- BuildListAdjacencyMatrices(PG, names, type=c('All','Shared', 'Specific'))
#' rowdata_Aggregation_sam(assay(Exp1_R25_pept,2), X.list)
rowdata_stats_Aggregation_sam <- function(qPepData, X.list){
  
  #X <- as.matrix(X)
  
  temp <- GetDetailedNbPeptidesUsed(X.list, qPepData)
  
  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", colnames(qPepData), sep="")
  rownames(pepSharedUsed) <- colnames(X.list$onlyShared)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", colnames(qPepData), sep="")
  rownames(pepSpecUsed) <- colnames(X.list$onlySpec)
  
  pepTotalUsed <- as.matrix(temp$nAll)
  colnames(pepTotalUsed) <- paste("pepTotal.used.", colnames(qPepData), sep="")
  rownames(pepTotalUsed) <- colnames(X.list$all)
  
  n <- GetDetailedNbPeptides(X.list)
  
  fd <- data.frame(colnames(X.list$all), 
                   nPepTotal = n$nTotal,
                   nPepShared = n$nShared, 
                   nPepSpec = n$nSpec, 
                   pepSpecUsed, 
                   pepSharedUsed, 
                   pepTotalUsed)


  # obj.prot@experimentData@other$typeOfData <-"protein"
  # #obj.prot <- addOriginOfValue(obj.prot)
  # obj.prot@experimentData@other$OriginOfValues <- NULL
  return (fd)
}




#' This function creates a column for the protein dataset after aggregation 
#' by using the previous peptide dataset.
#' 
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#' @param pepMetadata A data.frame of meta data of peptides. It is the fData 
#' of the MSnset object.
#' @param names The name of the column in fData(peptides_MSnset) that 
#' the user wants to keep in the new protein data.frame.
#' @param X The adjacency matrix used to agregate the peptides data.
#' @return A vector
#' @author Samuel Wieczorek
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' PG <- rowData(Exp1_R25_pept[['original_log']])[,metadata(Exp1_R25_pept)$parentProtId]
#' names <- names((Exp1_R25_pept[['original_log']]))
#' X.list <- BuildListAdjacencyMatrices(PG, names, type=c('All','Shared', 'Specific'))
#' X <- X.list$all
#' ft <- aggMetadata_sam(rowData(Exp1_R25_pept[['original_log']]), c('Sequence', 'Proteins'), X)
#' @export
aggMetadata_sam <- function(pepMetadata, names, X, simplify=TRUE){
  
  if(length(simplify)==1){
    if(isTRUE(simplify))
       simplify <- rep(TRUE,length(names))
    else 
      simplify <- rep(FALSE,length(names))
  } else if (length(simplify) != length(names)){
    warning('The length of parameter simplify must be equal to the length of the vector names. Abort.')
    return(NULL)
  }
  
  
  nbProt <- ncol(X)
  res <- setNames(data.frame(matrix(ncol = length(names), nrow = 0), stringsAsFactors = FALSE), names)
  proteinNames <- colnames(X)
  
  for (i in 1:length(proteinNames)){
    listeIndicePeptides <- names(which(X[,proteinNames[i]] == 1))
    listeData <- pepMetadata[listeIndicePeptides,names]
    listeData <- lapply(listeData, function(x){ gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(x, ',')))})
    
    line <- mapply(function(listeData, simplify) {if(isTRUE(simplify)) paste0(unique(listeData),collapse=', ') else paste0(listeData,collapse=', ')}, listeData, simplify)
    res <- setNames(rbind(res, data.frame(as.list(line))), names)
  }
  return(res)
}

##########################################################################################################################



#' This function creates a column for the protein dataset after agregation 
#' by using the previous peptide dataset. It is a parallel version of the function
#' \code{BuildColumnToProteinDataset}
#' 
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#' @param peptideData A data.frame of meta data of peptides. It is the fData 
#' of the MSnset object.
#' @param names The name of the column in fData(peptides_MSnset) that 
#' the user wants to keep in the new protein data.frame.
#' @param X The adjacency matrix used to agregate the peptides data.
#' @param simplify The names of the protein in the new dataset (i.e. rownames)
#' @return A data frame
#' @author Samuel Wieczorek
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:1000]
#' M <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' data <- Biobase::fData(obj.pep)
#' protData <- DAPAR::aggregateSum(obj.pep, M)
#' name <- "Protein_group_IDs"
#' proteinNames <- rownames(Biobase::fData(protData))
#' BuildColumnToProteinDataset_par(data, M, name,proteinNames )
#' }
#' @export
#' @import foreach
#' @importFrom doParallel registerDoParallel 
aggMetadata_parallel_sam <- function(peptideData, names, X, simplify=TRUE){
  doParallel::registerDoParallel()
  
  if(length(simplify)==1){
    if(isTRUE(simplify))
       simplify <- rep(TRUE,length(names))
       else 
         simplify <- rep(FALSE,length(names))
  } else if (length(simplify) != length(names)){
    warning('The length of parameter simplify must be equal to the length of the vector names. Abort.')
    return(NULL)
  }
  
  nbProt <- ncol(X)
  res <- setNames(data.frame(matrix(ncol = length(names), nrow = 0), stringsAsFactors = FALSE), names)
  proteinNames <- colnames(X)
  

  res <- foreach (i=1:length(proteinNames), .combine=rbind) %dopar% {
    listeIndicePeptides <- names(which(X[,proteinNames[i]] == 1))
    listeData <- pepMetadata[listeIndicePeptides,names]
    listeData <- lapply(listeData, function(x){ gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(x, ',')))})
    
    line <- mapply(function(listeData, simplify) {if(isTRUE(simplify)) paste0(unique(listeData),collapse=', ') else paste0(listeData,collapse=', ')}, listeData, simplify)
    setNames(rbind(res, data.frame(as.list(line))), names)
  }
  return(res)
}


