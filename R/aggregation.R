

#' @title Aggregate an assay's quantitative features which take into account
#' the peptides shared between proteins
#'
#' @description
#' 
#' This function aggregates the quantitative features of an assay,
#' applying an aggregation function (`fun`) to sets of features as
#' defined by the `fcol` feature variable. The new assay's features
#' will be named based on the unique `fcol` values.
#' This function is largely inspired by xxxx . The difference is that it can take into account the peptides shared between proteins.
#'
#'
#' @param object An instance of class [QFeatures].
#'
#' @param i The index or name of the assay which features will be
#'     aggregated the create the new assay.
#'
#'
#' @param aggType The type of peptides used for the aggregation. Possibla values are: 'all', 'onlyShared' and 'onlySPec'. This argument automatically
#' selects the corresponding adjacency matrix.
#' 
#' @param name A `character(1)` naming the new assay. Default is `newAssay`. Note that the function will fail if there's
#'     already an assay with `name`.
#'     
#' @param meta.names A vector of character strings that are the metadata of the peptides which needs to be aggregated
#' and kept in the protein dataset
#'
#' @param fun A function used for quantitative feature
#'     aggregation. See Details for examples.
#'
#' @param ... Additional parameters passed the `fun`.
#'
#' @return A `QFeatures` object with an additional assay.
#'
#' @details
#'
#' Aggregation is performed by a function that takes a matrix as
#' input and returns a xxxxx. Examples
#' thereof are
#'
#' - [Dapar2:aggSum()] to use the sum of each column (default);
#'
#' - [Dapar2:aggMean()] to use the sum of each column;
#'
#' - [Dapar2:aggIter()] to use the mean of each column;
#'
#' - [Dapar2:aggIterParallel()] same as previous function but use parallelism.
#'
#' - [Dapar2::aggTopn] to use the sum of each column;
#'
#' 
#' @seealso The *QFeatures* vignette provides an extended example and
#'     the *Processing* vignette, for a complete quantitative
#'     proteomics data processing pipeline.
#'
#' @name aggregateFeatures_sam
#'
#' @importFrom MsCoreUtils aggregate_by_vector robustSummary
#' @importFrom S4Vectors Hits
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' aggregateFeatures_sam(obj,2, aggType= 'all', name='aggregated', 
#' meta.names = 'Sequence', fun='aggTopn', n=3)
#' 
#' @export aggregateFeatures_sam
#' 
##-------------------------------------------------------------------------
# setMethod("aggregateFeatures_sam", "QFeatures",
#           function(object, i, X, name = "newAssay",
#                    fun = aggSum, ...)
#             .aggregateFeatures_sam(object, i, fcol, name, fun, ...))


aggregateFeatures_sam <- function(object, i, aggType='all', name, meta.names = NULL, fun, ...) {
  if (isEmpty(object))
    return(object)
  if (name %in% names(object))
    stop("There's already an assay named '", name, "'.")
  if (missing(aggType))
    stop("'aggType' is required.")    
  if (missing(i))
    stop("'i' is required.")  
  
  
  argg <- c(as.list(environment()), list(...))
  
  ## if the adjacency matrices are already present in the metadata of the SE, then load it
  ## else build it
  
  assay_i <- assay(object, i)
  object_i <- object[[i]]
  rowdata_i <- rowData(object[[i]])
  
  if (!(aggType %in% c('all', 'onlyShared', 'onlySpec'))){
    stop("Available values for 'typeMatAdj' are 'all', 'onlyShared', 'onlySpec'.")
  }
  
  if (is.null(metadata(object_i)$list.matAdj) || is.null(metadata(object_i)$list.matAdj[[aggType]])){
    ## by default, all three types of matrices are computed once for all and stored in the metadata slot of the SE. 
    ## This reduces the future computations
    object <- addListAdjacencyMatrices(object, i)
  }
  
  X <- metadata(object[[i]])$list.matAdj[[aggType]]
  
  
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
  
  aggregated_assay <- aggregate_with_matAdj(assay_i, X, fun, ...)
  
  
  # aggregated_rowdata <- Features::reduceDataFrame(rowdata_i, rowdata_i[[fcol]],
  #                                                 simplify = TRUE, drop = TRUE,
  #                                                 count = TRUE)
  
  # This part build general statistics about the aggregation
  aggregated_rowdata <- rowdata_stats_Aggregation_sam(assay_i, X)
  
  ## This part make the aggregation of the metadata of peptides which have been given
  ## by the parameter meta.names. It is column-binded to the previous aggregated_rowdata
  if (!is.null(meta.names)){
    aggregated_rowdata <- cbind(aggregated_rowdata, aggMetadata_sam(rowdata_i, meta.names, X))
  }
  
  se <- SummarizedExperiment(aggregated_assay,
                             rowData = aggregated_rowdata[rownames(aggregated_assay), ])
  metadata(se)$Params <- argg[-match(c('object', 'i', 'name', 'meta.names'), names(argg))]
  se <- SetTypeDataset(se, 'protein')
  
  # hits <- findMatches(rownames(aggregated_assay), groupBy)
  # rownames(aggregated_assay) : correspond au nom des proteines nouvellement creees
  # groupBy : correspond aaux id des peptides d'origine ?
  #hits <- findMatches(rownames(aggregated_assay), groupBy)
  ## from et to forment la definition du graphe pepetide-protein
  
  
  ## The following code is used to build the hits object originally built with findMatches
  ## in the class QFeatures. The vectors from and to are built explicitly with the 'which' function
  test <- which(as.matrix(X)==1, arr.ind=TRUE)
  from <- test[,'col']
  to <- test[,'row']
  hits <- S4Vectors::Hits(from=from, 
                          to=to, 
                          nLnode=ncol(X), 
                          nRnode=nrow(X),
                          sort.by.query=TRUE
  )
  
  elementMetadata(hits)$names_from <- rownames(assay_i)[hits@to]
  elementMetadata(hits)$names_to <- colnames(X)[hits@from]
  
  
  assayLinks <- AssayLink(name = name,
                          from = ifelse(is.character(i), i, names(object)[i]),
                          fcol = metadata(object)$parentProtId,
                          hits = hits)
  
  addAssay(object,
           se,
           name = name,
           assayLinks = assayLinks)
  
}




#' Method to create a list of three binary matrices with proteins in columns and peptides in lines.
#' 
#' @title Function matrix of appartenance group
#' 
#' @description All the matrices have the same dimensions. The first matrix correspond to all the peptides, the second one contains 
#' only the peptides shared with several proteins and the last matrix contains only the specific peptides. These matrices are useful 
#' is one aggregate with only a certain type of peptides or all of them. The list of matrices is stored in the slot 'xxx' of the metadata of the argument
#' obj (class 'SummarizedExperiment')
#' 
#' @param obj.se An object of class 'SummarizedExperiment' 
#' 
#' @return A list of three adjacency matrices
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' ll.x <- ComputeAdjacencyMatrices(obj[[2]])
#' 
#' @export
#' 
#' @importFrom Matrix Matrix
#' @import QFeatures
#' 
ComputeAdjacencyMatrices <- function(obj.se){
  if(class(obj.se) != 'SummarizedExperiment')
    stop("'obj.se' is not a 'SummarizedExperiment' object")

  md <- metadata(obj.se)
  if (is.null(md$parentProtId) || md$parentProtId == '' || nchar(md$parentProtId) == 0){
    warning("'parentProtId' is missing.")
    return(obj.se)
  } else if (!(md$parentProtId %in% colnames(rowData(obj.se)))){
    warning("'parentProtId' is not correctly set and does not seem to belongs to the dataset.")
    return(obj.se)
  }
  
  
  # A vector of proteins ids. The length of this vector is equal to the number of peptides one wants to aggregate, 
  # each line of it correspond to a peptide. Each element of this vector is either one element or a combination of elements seperated by a comma.
  plist <- rowData(obj.se)[,md$parentProtId]
  
  # A vector of names of peptides (a unique Id). The size of this vector is equal to the size of the parameter 'plist'.
  names <- names(obj.se)
  
  
  if (length(plist) != length(names(obj.se))){
    warning("The parameters 'plist' and 'names(obj.se)' must have the same length. Abort.")
    return(NULL)
  }
  
  Xshared <- Xspec <- Xall <- NULL
  
  # the separator that is allowed is ','
  pg <- gsub(";", ",", as.character(plist), fixed=TRUE)
  PG.l <- strsplit(as.character(pg), split=",", fixed=TRUE)
  t <- table(data.frame(A = rep(seq_along(PG.l), lengths(PG.l)),
                        B = unlist(PG.l)))
  
  #compute the matrix with all peptides
    Xall <- Matrix::Matrix(t, sparse=T, dimnames = list(names(obj.se), colnames(t)))

  #compute the matrix with only shared peptides
    tmpShared <- t
    ll <- which(rowSums(tmpShared)==1)
    if (length(ll) > 0)
      tmpShared[ll,] <- 0

    Xshared <- Matrix::Matrix(tmpShared, sparse=T, dimnames = list(names(obj.se), colnames(t)))

  #compute the matrix with only specific peptides
    tmpSpec <- t
    ll <- which(rowSums(tmpSpec)>1)
    if (length(ll) > 0)
      tmpSpec[ll,] <- 0

    Xspec <- Matrix::Matrix(tmpSpec, sparse=T, dimnames = list(names(obj.se), colnames(t)))

    
  x.list <- list(all = Xall, 
                 onlyShared = Xshared, 
                 onlySpec = Xspec)
  
  
  return (x.list)
}




#' This function computes few values about the adjacency matrix such as the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two. 
#' 
#' @title Computes the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two.
#' 
#' @param X The adjacency matrix with both specific and shared peptides.
#' 
#' @return A list of values:
#' * nbPeptides: the number of peptides in the matrix,
#' nbSpecificPeptides: the number of specific peptides in the matrix,
#' nbSharedPeptides: the number of shared peptides in the matrix,
#' nbProt: the number of proteins in the matrix,
#' protOnlyUniquePep: the list of proteins only defined by specific peptides,
#' protOnlySharedPep: the list of proteins only defined by shared peptides,
#' protMixPep: the list of proteins defined by both shared and specific peptides.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' matAdjStats(X)
#' 
#' @export
#' 
matAdjStats <- function(X){
  if (is.null(X)){
    warning('The adjacency matrix is NULL.')
    return(NULL)
  }
  
  
  ind.shared.Pep <- which(rowSums(as.matrix(X))>1)
  ind.unique.Pep <- which(rowSums(as.matrix(X))==1)
  
  M.shared.Pep <- X[ind.shared.Pep,]
  M.shared.Pep <- M.shared.Pep[,-which(colSums(as.matrix(M.shared.Pep))==0)]
  
  M.unique.Pep <- X[ind.unique.Pep,]
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' GraphPepProt_hc(X)
#' 
#' @import highcharter
#' 
#' @export
#' 
GraphPepProt_hc <- function(X, type = 'all'){
  if (is.null(X)){
    warning('The parameter mat is empty.')
    return (NULL)
  } 

  t <- t(X)
  t <- apply(X, 2, sum, na.rm=TRUE)
  tab <- table(t)
  conds <- names(tab)
  
  h1 <-  highchart() %>%
    dapar_hc_chart(chartType = "column") %>%
    hc_title(text = paste0("Distribution of ",type, " peptides w.r.t. proteins")) %>%
    hc_add_series(data=tab,type="column", colorByPoint = TRUE) %>%
    hc_colors('orange') %>%
    hc_plotOptions( column = list(stacking = "normal"),
                    animation=list(duration = 100)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = conds, title = list(text = "Number of peptides")) %>%
    hc_yAxis(categories = conds, title = list(text = "Number of proteins")) %>%
    dapar_hc_ExportMenu(filename = "HistoMatAdj") %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "{point.y}")
  
  h1
  
}




#' Method to compute the number of quantified peptides used for aggregating each protein
#' 
#' @title Computes the number of peptides used for aggregating each protein
#' 
#' @param qPepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @return A data.frame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' qPepData <- assay(obj,2)
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' n <- GetNbPeptidesUsed(qPepData, X)
#' 
#' @export
#' 
GetNbPeptidesUsed <- function(qPepData, X){

  qPepData[!is.na(qPepData)] <- 1
  qPepData[is.na(qPepData)] <- 0
  
  pep <- t(X) %*% qPepData
  
  return(pep)
}




#' Method to compute the detailed number of quantified peptides used for aggregating each protein
#' 
#' @title Computes the detailed number of peptides used for aggregating each protein w.r.t NA values in 
#' peptide quantitative dataset. Even if a peptide is part of a protein, if its value is NA, it do not be used
#' for aggreation.
#' 
#' @param X An adjacency matrix.
#' 
#' @param qPepData A data.frame of quantitative data
#' 
#' @seealso The function 'addListAdjacencyMatrices'.
#' 
#' @return A list of three items:
#' * nAll: the number of pall eptides used for aggregation,
#' * nShared: the number of shared peptides used for aggregation,
#' * nSpec: the total number of pspecific eptides used for aggregation.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' qPepData <- assay(obj,2)
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' n <- GetDetailedNbPeptidesUsed(X, qPepData)
#'  
#' @export
#' 
GetDetailedNbPeptidesUsed <- function(X, qPepData){
  
  res <- NULL
  
  qPepData[!is.na(qPepData)] <- 1
  qPepData[is.na(qPepData)] <- 0

  res <- t(as.matrix(X)) %*% qPepData
  
  return(res)
}


#' Method to compute the detailed number of peptides for each protein
#' 
#' @title Computes the detailed number of peptides for each protein 
#' 
#' @param X An adjacency matrices
#' 
#' @return A data.frame containing the following vectors:
#' * nTotal: The number of peptides for each protein,
#' * nShared: the number of shared peptides for each protein,
#' * nSPec: the number of specific peptides for each protein.
#' 
#' Each row of these three vectors represent a protein. Thus, the length of these vectors is equal 
#' to the number of proteins.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' n <- GetDetailedNbPeptides(X)
#' 
#' @export
#' 
GetDetailedNbPeptides <- function(X){
  
  n <- rowSums(t(as.matrix(X)))
  
  return(n)
  
}





#' Method to aggregate peptides into proteins with the sum of the quantitative data per conditions.
#' 
#' @title aggregate peptides into proteins with the sum of the quantitative data per conditions.
#'  
#' @param qPepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @return A matrix containing the quantitative aggregated data for proteins.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' n <- inner.sum(assay(obj[[2]]), X)
#' 
#' @export
#' 

inner.sum <- function(qPepData, X){
  qPepData[is.na(qPepData)] <- 0
  
  Mp <- t(X) %*% qPepData
  return(Mp)
}



#' Method to aggregate peptides into proteins with the mean of the quantitative data per conditions.
#' 
#' @title Aggregate peptides into proteins with the mean of the quantitative data per conditions. 
#' 
#' @param qPepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @return A matrix containing the quantitative aggregated data for proteins.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' inner.mean(assay(obj[[2]]), X)
#' 
#' @export
#' 
inner.mean <- function(qPepData, X){
  
  X <- as.matrix(X)
  
  Mp <- inner.sum(qPepData, X)
  Mp <- Mp / GetNbPeptidesUsed(qPepData, X)
  
  return(Mp)
}





#' Method to aggregate peptides into proteins with the top n approach.
#' 
#' @title Method to aggregate peptides into proteins with the top n approach.
#' 
#' @param qPepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @param method Method used to aggregate (see function xxx)
#' 
#' @param n An integer which is the number of top peptides used for aggregation.
#' 
#' @return A matrix containing the quantitative aggregated data for proteins.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' inner.aggregate.topn(assay(obj[[2]]), X, n=3)
#' 
#' @importFrom stats median
#' 
#' @export
#' 
inner.aggregate.topn <-function(qPepData, X, method='Mean', n=10){
  
  X <- as.matrix(X)
  
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' qPepData <- assay(obj[[2]])
#' inner.aggregate.iter(qPepData, X)
#' 
#' @export
#' 
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
  }
  return(as.matrix(yprot))
}




#' This function computes the intensity of proteins based on the sum of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins with the sum of the intensities
#' of their peptides.
#' 
#' @param qPepData A matrix of intensities of peptides.
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A matrix of intensities of proteins.
#' 
#' @author Alexia Dorffer, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]], 'all'))
#' aggSum(assay(obj[[2]]), X)
#' 
#' @export
#' 
aggSum <- function(qPepData, X){
  #qPepData <- 2^(qPepData)
  protData <- inner.sum(qPepData, X)

  return(protData)
}



#' This function computes the intensity of proteins as the mean of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins as the mean of the intensities
#' of their peptides.
#' 
#' @param qPepData A matrix of intensities of peptides.
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Alexia Dorffer
#' 
#' @examples 
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' aggMean(assay(obj[[2]]), X)
#' 
#' @export
#' 
aggMean <- function(qPepData, X){
  #qPpepData <- 2^(qPepData)
  protData <- inner.mean(qPepData, X)
  return(protData)
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
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
aggIterParallel <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
  if (is.null(conditions)){
    warning('The parameter \'conditions\' is NULL: the aggregation cannot be process.')
    return(NULL)
  }
  doParallel::registerDoParallel()
  
  #qPepData <- 2^(qPepData)
  protData <- matrix(rep(0,ncol(X)*nrow(X)), nrow=ncol(X))
  cond <- NULL
  protData <- foreach(cond = 1:length(unique(conditions)),
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' conditions <- colData(obj)$Condition
#' aggIter(assay(obj,2), X, conditions)
#' 
#' @export
aggIter <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
  if (is.null(conditions)){
    warning("The parameter 'conditions' is NULL: the aggregation cannot be process.")
    return(NULL)
  }
  
  #qPepData <- 2^(qPepData)
  
  protData <- matrix(rep(0,ncol(X)*ncol(qPepData)), nrow=ncol(X))
  
  for (cond in unique(conditions)){
    condsIndices <- which(conditions == cond)
    qData <- qPepData[,condsIndices]
    protData[,condsIndices]  <- inner.aggregate.iter(qData, X, init.method, method, n)
  }
  
  return(protData)

}




#' This function computes the intensity of proteins as the sum of the 
#' intensities of their n best peptides.
#' 
#' @title Compute the intensity of proteins as the sum of the intensities of their n best peptides.
#' 
#' @param qPepData A matrix of intensities of peptides
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @param method Default is 'Mean'
#' 
#' @param n The maximum number of peptides used to aggregate a protein.
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Alexia Dorffer, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' aggTopn(assay(obj,2), X, n=3)
#' 
#' @export
#' 
aggTopn <- function(qPepData, X,  method='Mean', n=10){
  #qPepData <- 2^(qPepData)
  
  protData <- inner.aggregate.topn(qPepData, X, method=method, n)
  return(protData)
}



#' @title Aggreagate quantitative features.
#' 
#' @description
#' This function takes a matrix of quantitative features `x` and a
#' factor (of length equal to `nrow(x)`) defining subsets, and
#' applies a user-defined function to aggregate each subset into a
#' vector of quantitative values. This function is the same as 'aggregate_by_vector' from the package QFeatures
#' and it is used with DAPAR aggregation functions
#'
#' User-defined functions must thus return a matrix of dimensions equal
#' to `X`. Examples thereof are
#'
##' - [Dapar2:aggSum()] to use the sum of each column (default);
##'
##' - [Dapar2:aggMean()] to use the sum of each column;
##'
##' - [Dapar2:aggIter()] to use the mean of each column;
##'
##' - [Dapar2:aggIterParallel()] same as previous function but use parallelism.
##'
##' - [Dapar2::aggTopn] to use the sum of each column;
#'
#' @param qPepData xxxx. 
#' 
#' @param X An adjacency matrix.
#' 
#' @param FUN A `function` to be applied to the subsets of `x`.
#' 
#' @param ... Additional arguments passed to `FUN`.
#' 
#' @return A new `matrix` of dimensions xxxx
#' 
#' @author Samuel Wieczorek
#'
#' @family Quantitative feature aggregation
#' 
#' @export
#' 
aggregate_with_matAdj <- function(qPepData, X, FUN, ...){
  res <- do.call(FUN, list(qPepData, X, ...))
  res
}








#' Method to finalize the aggregation process
#' 
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' rowdata_stats_Aggregation_sam(assay(obj,2), X)
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




#' This function aggregates features of metadata of peptides into a dataframe that is to be included
#' in the rowdata of the new protein dataset.
#' 
#' @title This function aggregates features of metadata of peptides into a dataframe that is to be included
#' in the rowdata of the new protein dataset.
#' 
#' @param pepMetadata A data.frame of meta data of peptides. It is the rowData 
#' of the QFeatures object.
#' 
#' @param names A vector of column names of the rowdata of peptides.
#' 
#' @param X The adjacency matrix used to agregate the peptides data.
#' 
#' @param simplify This parameter is either a boolean or a vector of boolean. The TRUE value means that, for a given protein,
#' duplicates metadata are reduced into one. The FALSE value means that, for a given protein, duplicates metadata are not
#' reduced: one keep the real original values.
#' If the parameter is a vector, its size must be equal to the size of the parameter names. Each item in 'simplify' is the indication
#' for the item of the same indice in the vector 'names' to apply the reduction or not.
#' 
#' @return A data frame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]], 'all'))
#' ft <- aggMetadata_sam(rowData(obj[[2]]), c('Sequence'), X)
#' 
#' @importFrom stats setNames
#' 
#' @export
#' 
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
  res <- stats::setNames(data.frame(matrix(ncol = length(names), nrow = 0), stringsAsFactors = FALSE), names)
  proteinNames <- colnames(X)
  
  for (i in 1:length(proteinNames)){
    listeIndicePeptides <- names(which(X[,proteinNames[i]] == 1))
    listeData <- pepMetadata[listeIndicePeptides,names]
    
    if (length(names)==1){
     # listeData <- setNames(DataFrame(matrix(pepMetadata[listeIndicePeptides,names])), names)
      listeData <-gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(listeData, ',')))
      if(isTRUE(simplify)) {
        line <- setNames(paste0(unique(listeData),collapse=', ') , names)
      } else {
        line <- setNames(paste0(listeData,collapse=', '), names)
      }
    } else {
    listeData <- pepMetadata[listeIndicePeptides,names]
    listeData <- lapply(listeData, function(x){ gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(x, ',')))})
    line <- mapply(function(listeData, simplify) {if(isTRUE(simplify)) 
                                                    paste0(unique(listeData),collapse=', ') 
                                                  else 
                                                    paste0(listeData,collapse=', ')}, listeData, simplify)
    
    }
    
    
    res <- stats::setNames(rbind(res, data.frame(as.list(line))), names)
  }
  return(res)
}




#' This function aggregates features of metadata of peptides into a dataframe that is to be included
#' in the rowdata of the new protein dataset. It is the same function as aggMetadata_sam but with parallelism.
#' 
#' @title This function aggregates features of metadata of peptides into a dataframe that is to be included
#' in the rowdata of the new protein dataset.
#' 
#' @param pepMetadata A data.frame of meta data of peptides. It is the rowData 
#' of the QFeatures object.
#' 
#' @param names A vector of column names of the rowdata of peptides.
#' 
#' @param X The adjacency matrix used to agregate the peptides data.
#' 
#' @param simplify This parameter is either a boolean or a vector of boolean. The TRUE value means that, for a given protein,
#' duplicates metadata are reduced into one. The FALSE value means that, for a given protein, duplicates metadata are not
#' reduced: one keep the real original values.
#' If the parameter is a vector, its size must be equal to the size of the parameter names. Each item in 'simplify' is the indication
#' for the item of the same indice in the vector 'names' to apply the reduction or not.
#' 
#' @return A data frame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' ft <- aggMetadata_parallel_sam(rowData(obj[[2]]), c('Sequence', 'Proteins'), X)
#' 
#' @export
#' 
#' @import foreach
#' @importFrom doParallel registerDoParallel 
#' @importFrom stats setNames
#' 
aggMetadata_parallel_sam <- function(pepMetadata, names, X, simplify=TRUE){
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
  res <- stats::setNames(data.frame(matrix(ncol = length(names), nrow = 0), stringsAsFactors = FALSE), names)
  proteinNames <- colnames(X)
  
  i <- NULL
  res <- foreach (i=1:length(proteinNames), .combine=rbind, .packages = c("S4Vectors","Matrix")) %dopar% {
    listeIndicePeptides <- names(which(X[,proteinNames[i]] == 1))
    listeData <- pepMetadata[listeIndicePeptides,names]
    
    if (length(names)==1){
      # listeData <- setNames(DataFrame(matrix(pepMetadata[listeIndicePeptides,names])), names)
      listeData <-gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(listeData, ',')))
      if(isTRUE(simplify)) {
        line <- setNames(paste0(unique(listeData),collapse=', ') , names)
      } else {
        line <- setNames(paste0(listeData,collapse=', '), names)
      }
    } else {
      listeData <- pepMetadata[listeIndicePeptides,names]
      listeData <- lapply(listeData, function(x){ gsub(pattern = "\\s", replacement = "", x = unlist(strsplit(x, ',')))})
      line <- mapply(function(listeData, simplify) {if(isTRUE(simplify)) 
        paste0(unique(listeData),collapse=', ') 
        else 
          paste0(listeData,collapse=', ')}, listeData, simplify)
    }
    
    stats::setNames(rbind(res, data.frame(as.list(line))), names)
  }
  return(res)
}

