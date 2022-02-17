
#' @title Aggregate an assay's quantitative features
#'
#' @description
#'
#' This function aggregates the quantitative features of an assay,
#' applying a summarisation function (`fun`) to sets of features.
#' The `fcol` variable name points to a rowData column that defines
#' how to group the features during aggregate. This variable can
#' eigher be a vector (we then refer to an *aggregation by vector*)
#' or an adjacency matrix (*aggregation by matrix*).
#'
#' @param object An instance of class [QFeatures] or [SummarizedExperiment].
#'
#' @param i The index or name of the assay which features will be
#'     aggregated the create the new assay.
#'
#' @param fcol A `character(1)` naming a rowdata variable (of assay
#'     `i` in case of a `QFeatures`) defining how to aggregate the
#'     features of the assay. This variable is either a `character`
#'     or a (possibly sparse) matrix. See below for details.
#'
#' @param name A `character(1)` naming the new assay. Default is
#'     `newAssay`. Note that the function will fail if there's
#'     already an assay with `name`.
#'
#' @param fun A function used for quantitative feature
#'     aggregation. See Details for examples.
#'
#' @param ... Additional parameters passed the `fun`.
#'
#' @return A `QFeatures` object with an additional assay or a
#'  `SummarizedExperiment` object (or subclass thereof).
#'
#' @seealso The *QFeatures* vignette provides an extended example and
#'     the *Processing* vignette, for a complete quantitative
#'     proteomics data processing pipeline. The
#'     [MsCoreUtils::aggregate_by_vector()] manual page provides
#'     further details.
#'
#' @aliases aggregateFeatures aggregateFeatures,QFeatures-method
#'     aggcounts aggcounts,SummarizedExperiment-method
#'     adjacencyMatrix,SummarizedExperiment-method
#'     adjacencyMatrix,QFeatures-method
#'
#' @name aggregateQmetadata
#'
#' @rdname DaparToolshed-aggregate
#'
#' @examples
#'
NULL

#' @exportMethod aggregateQmetadata
#' @rdname qMetadata-aggregate
setMethod("aggregateQmetadata", "QFeatures",
          function(object, i, fcol, name = "newAssay",
                   fun = aggQmeta, ...) {
            if (isEmpty(object))
              return(object)
            if (!(name %in% names(object)))
              stop("An assay named '", name, "' is not found.")
            if (missing(i))
              i <- main_assay(object)
            
            
            # Aggregate the quantitative metdata
            aggQ <- aggregateQmetadata(object[[i]], fun, conds = colData(object)$Condition)
            
            # Add the aggregated qMetadata to the se
            qMetadata(object[[name]]) <- aggQ
            
            object
          })


#' @exportMethod aggregateQmetadata
#' @rdname qMetadata-aggregate
#' @export
setMethod("aggregateQmetadata", "SummarizedExperiment",
          function(object, fun , conds)
            do.call(fun, list(object,  conds))
          )









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


#' @title Finalizes the aggregation process 
#' 
#' @param qPepData A data.frame of quantitative data not logged of peptides
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A protein object of class \code{SummarizedExperiment}
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @examples
#' library(QFeatures)
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' rowdata_stats_Aggregation_sam(assay(obj,2), X)
#' 
#' @rdname DaparToolshed-aggregate
#' 
rowdata_stats_Aggregation_sam <- function(qPepData, X){
  
  X <- as.matrix(X)
  
  temp <- GetDetailedNbPeptidesUsed(X, qPepData)
  
  pepUsed <- as.matrix(temp)
  colnames(pepUsed) <- paste("pep_used_", colnames(qPepData), sep="")
  rownames(pepUsed) <- colnames(X)
  
  n <- GetDetailedNbPeptides(X)
  
  fd <- data.frame(colnames(X), 
                   n, 
                   pepUsed)
  
  return (fd)
}







#' @title Get the type of dataset
#' @description xxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[1:100]
#' X <- adjacencyMatrix(obj[[2]])
#' agg.qMmetadata <- aggregateQmetadata(obj, 2, X)
#' agg.meta <- aggregateQmetadata(obj[[2]], X)
#' 
#' @export
#' @return NA
#' 
NULL

#' 
#' @param object  An object of class 'SummarizedExperiment'
#' @param conds xxx
#' 
#' @return NA
#' 
#' @rdname DaparToolshed-aggregate
#' 
aggQmeta <- function(object, conds) {
  stopifnot(inherits(object, "SummarizedExperiment"))
  
  qMeta = qMetadata(object)
  level = typeDataset(object)
  X <- adjacencyMatrix(object)
  
  rowcol <- function(meta.col, X.col)
    meta.col[X.col > 0]
  
  df <- data.frame(stringsAsFactors = TRUE)
  for (j in 1:ncol(qMeta))
    for(i in 1:ncol(X))
      df[i, j] <- qMetadata_combine( rowcol(qMeta[,j], X[,i]), 
                                     level)
  
  df[df=='NA'] <- NA
  dimnames(df) <- list(colnames(X), colnames(qMeta))
  # Delete protein with only NA
  
  # Post processing of metacell to discover 'imputed POV', 'imputed MEC'
  df <- Set_POV_MEC_tags(conds, df, level)
  
  df
}




#' @title Finalizes the aggregation process 
#' 
#' @param obj.pep A peptide object of class \code{MSnset}
#' 
#' @param pepData xxxx
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @param protData xxxxx
#' 
#' @param protMetacell xxx
#' 
#' @return A protein object of class \code{MSnset}
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @importFrom utils installed.packages
#' 
#' @examples
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[1:100]
#' FinalizeAggregation(obj, 2)
#' 
#' @rdname DaparToolshed-aggregate
#' 
setMethod("FinalizeAggregation", "QFeatures",
          function(object, from, to, ...){
            
            if (isEmpty(object))
              return(object)

            if (missing(from))
              stop("Provide index or name of assay to be processed")
            if (length(from) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(from)) from <- names(object)[[from]]
            
            if (missing(to))
              stop("Provide index or name of assay to be processed")
            if (length(to) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(to)) to <- names(object)[[to]]
            
            
            
  
  object.from <- object[[from]]
  object.to <- object[[to]]
  
  from.data <- assay(object.from)
  to.data <- assay(object.to)
  
  to.data <- as.matrix(to.data)
  to.data[to.data==0] <- NA
  to.data[is.nan(to.data)] <- NA
  to.data[is.infinite(to.data)] <-NA
  
  X.all <- adjacencyMatrix(object[[from]])
  X.spec <- updateAdjacencyMatrix(X.all, mode = 'onlySpec')
  X.shared <- updateAdjacencyMatrix(X.all, mode = 'onlyShared')
  
  rowData(object[[to]])$pepSpecUsed = t(as.matrix(X.spec)) %*% (!is.na(qData))
  rowData(object[[to]])$pepSharedUsed = t(as.matrix(X.shared)) %*% (!is.na(qData))
  rowData(object[[to]])$pepTotalUsed = t(as.matrix(X.all)) %*% (!is.na(qData))

  
  metadata(object[[to]]) <- metadata(object[[from]])
  typeDataset(object[[to]]) <- "protein"
  idcol(object[[to]]) <- 'proteinId'
  metadata(object[[to]])$proteinId <- 'proteinId'
   
  tryCatch({
    find.package("Prostar")
    find.package("DaparToolshed")
    
    version(object[[to]]) <- list(Prostar = Biobase::package.version('Prostar'),
                                   DaparToolshed_Version = Biobase::package.version('DaparToolshed')
    )
  },
  error = function(e) {
    version(object[[to]]) <- list(Prostar_Version = NA,
                                  DaparToolshed_Version = NA)
  }
  )
  
  return (object)
          }
)

