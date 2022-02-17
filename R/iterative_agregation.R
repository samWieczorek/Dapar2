




#' Method to aggregate peptides into proteins with the iterative approach.
#' 
#' @title Method to aggregate peptides into proteins with the iterative approach.
#'  
#' @param qPepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' 
#' @param method The method used for xxx. Default value is 'Mean'.
#' 
#' @param n An integer which is xxx
#' 
#' @return A  matrix containing the quantitative aggregated data for proteins.
#' 
#' @author Samuel Wieczorek, Thomas Burger
#' 
#' @examples
#' library(QFeatures)
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' qPepData <- assay(obj[[2]])
#' inner.aggregate.iter(qPepData, X)
#' 
#' @rdname DaparToolshed-aggregate
#' 
inner.aggregate.iter <- function(qData, 
                                 X, 
                                 init.method = 'colSumsMat', 
                                 iter.method = 'Mean', 
                                 n = NULL){
  
  if (!(init.method %in% c("colSumsMat", "colMeansMat"))) {
    warning("Wrong parameter init.method")
    return(NULL)
  }
  
  if (!(iter.method %in% c("onlyN", "Mean"))){
    warning("Wrong parameter iter.method")
    return(NULL)
  }
  
  
  if (iter.method=='onlyN' && is.null(n)){
    warning("Parameter n is null")
    return(NULL)
  }
  
  yprot <- do.call(init.method, list(X, qData))
  
  conv <- 1
  
  while(conv > 10**(-10)){
    mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot[is.na(mean.prot)] <- 0
    
    X.tmp <- mean.prot*X
    X.new <- X.tmp/rowSums(as.matrix(X.tmp), na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    # l'appel a la fonction ci-dessous depend des parametres choisis par l'utilisateur
    switch(iter.method,
           Mean = yprot <- colMeansMat(X.new, qData),
           onlyN = yprot <- inner.aggregate.topn(qData, X.new, 'Mean', n)
    )
    
    mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot.new[is.na(mean.prot.new)] <- 0
    
    conv <- mean(abs(mean.prot.new - mean.prot))
  }
  return(as.matrix(yprot))
}



#' Method to aggregate peptides into proteins with the iterative approach with use of parallelism.
#' 
#' @title Method to aggregate peptides into proteins with the iterative approach with use of parallelism.
#'  
#' @param qPepData A matrix of intensities of peptides.
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @param conditions xxxx
#' 
#' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' 
#' @param method The method used for xxx. Default value is 'Mean'.
#' 
#' @param n xxxx
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' library(doParallel)
#' library(QFeatures)
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' conditions <- SummarizedExperiment::colData(obj)$Condition
#' aggIterParallel(assay(obj,2), X, conditions)
#' 
#' @export
#' 
#' @import doParallel
#' @import foreach
#' 
#' @rdname DaparToolshed-aggregate
#' 
aggIterParallel <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
  if (is.null(conditions)){
    warning('The parameter \'conditions\' is NULL: the aggregation cannot be process.')
    return(NULL)
  }
  doParallel::registerDoParallel()
  
  #qPepData <- 2^(qPepData)
  protData <- matrix(rep(0,ncol(X)*nrow(X)), nrow=ncol(X))
  cond <- NULL
  protData <- foreach(cond = seq_len(length(unique(conditions))),
                      .combine=cbind,
                      .export=c("inner.aggregate.iter", "inner.sum", "inner.mean","GetNbPeptidesUsed"),
                      .packages = "Matrix") %dopar% {
                        condsIndices <- which(conditions == unique(conditions)[cond])
                        qData <- qPepData[,condsIndices]
                        inner.aggregate.iter(qData, X, init.method, method, n) 
                      }
  
  protData <- protData[,colnames(qPepData)]
  return(protData)
}




#' Method to aggregate peptides into proteins with the iterative approach 
#' 
#' @title Method to aggregate peptides into proteins with the iterative approach.
#'  
#' @param qPepData A matrix of intensities of peptides.
#' 
#' @param X An adjacency matrix
#' 
#' @param conditions xxx
#' 
#' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' 
#' @param method The method used for xxx. Default value is 'Mean'.
#' 
#' @param n xxxx
#' 
#' @return A matrix of protein intensities.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' conditions <- colData(obj)$Condition
#' aggIter(assay(obj,2), X, conditions)
#' 
#' @export
#' 
#' @rdname DaparToolshed-aggregate
#'
aggIterative <- function(x,
                         MAT, 
                         conditions = NULL, 
                         init.method = 'Sum', 
                         iter.method = 'Mean', 
                         n = NULL){
  if (is.null(conditions)){
    warning("The parameter 'conditions' is NULL: the aggregation cannot be process.")
    return(NULL)
  }
  
  protData <- matrix(rep(0,ncol(x)*ncol(MAT)), nrow=ncol(x))
  
  for (cond in unique(conditions)){
    ind <- which(conditions == cond)
    protData[,ind]  <- inner.aggregate.iter(MAT[,ind],
                                            x, 
                                            init.method,
                                            iter.method, 
                                            n)
  }
  
  return(protData)
  
}

