
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
#' - [DaparToolshed:aggSum()] to use the sum of each column (default);
#'
#' - [DaparToolshed:aggMean()] to use the sum of each column;
#'
#' - [DaparToolshed:aggIter()] to use the mean of each column;
#'
#' - [DaparToolshed:aggIterParallel()] same as previous function but use parallelism.
#'
#' - [DaparToolshed::aggTopn] to use the sum of each column;
#'
#' 
#' @seealso The *QFeatures* vignette provides an extended example and
#'     the *Processing* vignette, for a complete quantitative
#'     proteomics data processing pipeline.
#'
#'
#' @importFrom MsCoreUtils aggregate_by_vector robustSummary
#' @importFrom S4Vectors Hits
#' @importFrom MultiAssayExperiment isEmpty
#'
#' @examples
#' feat2 <- readRDS('~/GitHub/DaparToolshedData/data/Exp2_R100_pept.rda')
#' feat2 <- feat2[1:10,]
#' # Builds the adjacency matrix w.r.t. 'mode'
#' X <- makeAdjacencyMatrix(rowData(feat2[['original_log']])[,'Protein_group_IDs'])
#' rownames(X) <- rownames(feat2[['original_log']])
#' rowData(feat2[['original_log']])[['adjacencyMatrix']] <- NULL
#' adjacencyMatrix(feat2[['original_log']]) <- CustomAdjMat(feat2[['original_log']], X, mode = 'all')
#' feat2 <- AggregateFeatures4Prostar(object = feat2, 
#'                                    i = 2, 
#'                                    fcol = "adjacencyMatrix", 
#'                                    fun = aggSum)
#' 
#' @export
#' 
setMethod("aggregateFeatures4Prostar", "QFeatures",
           function(object, i, name = "newAssay", fcol, fun, ...){
  object <- aggregateFeatures(object,
                              i, 
                              name = name, 
                              fcol = "adjacencyMatrix", 
                              fun = fun,
                              ...)
  
  # Aggregate quantitative cell metadata
  object <- aggregateQmetadata(object,
                               i, 
                               name = name, 
                               fcol = "qMetadata",
                               fun = aggQmeta)
  
  object
}
)







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
#' The rowData of the aggregated `SummarizedExperiment` assay
#' contains a `.n` variable that provides the number of parent
#' features that were aggregated.
#'
#' When aggregating with a vector, the newly aggregated
#' `SummarizedExperiment` assay also contains a new `aggcounts` assay
#' containing the aggregation counts matrix, i.e. the number of
#' features that were aggregated for each sample, which can be
#' accessed with the `aggcounts()` accessor.
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
#' @details
#'
#' Aggregation is performed by a function that takes a matrix as
#' input and returns a vector of length equal to `ncol(x)`. Examples
#' thereof are
#'
#' - [MsCoreUtils::medianPolish()] to fits an additive model (two way
#'   decomposition) using Tukey's median polish_ procedure using
#'   [stats::medpolish()];
#'
#' - [MsCoreUtils::robustSummary()] to calculate a robust aggregation
#'   using [MASS::rlm()] (default);
#'
#' - [base::colMeans()] to use the mean of each column;
#'
#' - `colMeansMat(x, MAT)` to aggregate feature by the calculating
#'    the mean of peptide intensities via an adjacency matrix. Shared
#'    peptides are re-used multiple times.
#'
#' - [matrixStats::colMedians()] to use the median of each column.
#'
#' - [base::colSums()] to use the sum of each column;
#'
#' - `colSumsMat(x, MAT)` to aggregate feature by the summing the
#'    peptide intensities for each protein via an adjacency
#'    matrix. Shared peptides are re-used multiple times.
#'
#' See [MsCoreUtils::aggregate_by_vector()] for more aggregation functions.
#'
#' @section Missing quantitative values:
#'
#' Missing quantitative values have different effects based on the
#' aggregation method employed:
#'
#' - The aggregation functions should be able to deal with missing
#'   values by either ignoring or propagating them. This is often
#'   done with an `na.rm` argument, that can be passed with
#'   `...`. For example, `rowSums`, `rowMeans`, `rowMedians`,
#'   ... will ignore `NA` values with `na.rm = TRUE`, as illustrated
#'   below.
#'
#' - Missing values will result in an error when using `medpolish`,
#'   unless `na.rm = TRUE` is used. Note that this option relies on
#'   implicit assumptions and/or performes an implicit imputation:
#'   when summing, the values are implicitly imputed by 0, assuming
#'   that the `NA` represent a trully absent features; when
#'   averaging, the assumption is that the `NA` represented a
#'   genuinely missing value.
#'
#' - When using robust summarisation, individual missing values are
#'   excluded prior to fitting the linear model by robust
#'   regression. To remove all values in the feature containing the
#'   missing values, use [filterNA()].
#'
#' More generally, missing values often need dedicated handling such
#' as filtering (see [filterNA()]) or imputation (see [impute()]).
#'
#' @section Missing values in the row data:
#'
#' Missing values in the row data of an assay will also impact the
#' resulting (aggregated) assay row data, as illustrated in the
#' example below. Any feature variables (a column in the row data)
#' containing `NA` values will be dropped from the aggregated row
#' data. The reasons underlying this drop are detailed in the
#' `reduceDataFrame()` manual page: only invariant aggregated rows,
#' i.e. rows resulting from the aggregation from identical variables,
#' are preserved during aggregations.
#'
#' The situation illustrated below should however only happen in rare
#' cases and should often be imputable using the value of the other
#' aggregation rows before aggregation to preserve the invariant
#' nature of that column. In cases where an `NA` is present in an
#' otherwise variant column, the column would be dropped anyway.
#'
#' @section Using an adjacency matrix:
#'
#' When considering non-unique peptides explicitly, i.e. peptides
#' that map to multiple proteins rather than as a protein group, it
#' is convenient to encode this ambiguity explicitly using a
#' peptide-by-proteins (sparse) adjacency matrix. This matrix is
#' typically stored in the rowdata and set/retrieved with the
#' [adjacencyMatrix()] function. It can be created manually (as
#' illustrated below) or using `PSMatch::makeAdjacencyMatrix()`.
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
#' @rdname QFeatures-aggregate
#'
#' @importFrom MsCoreUtils aggregate_by_vector aggregate_by_matrix robustSummary colCounts
#'
#' @examples
#'
#' ## ---------------------------------------
#' ## An example QFeatures with PSM-level data
#' ## ---------------------------------------
#' data(feat1)
#' feat1
#'
#' ## Aggregate PSMs into peptides
#' feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
#' feat1
#'
#' ## Aggregate peptides into proteins
#' feat1 <- aggregateFeatures(feat1, "peptides", "Protein", name = "proteins")
#' feat1
#'
#' assay(feat1[[1]])
#' assay(feat1[[2]])
#' aggcounts(feat1[[2]])
#' assay(feat1[[3]])
#' aggcounts(feat1[[3]])
#'
#' ## --------------------------------------------
#' ## Aggregation with missing quantitative values
#' ## --------------------------------------------
#' data(ft_na)
#' ft_na
#'
#' assay(ft_na[[1]])
#' rowData(ft_na[[1]])
#'
#' ## By default, missing values are propagated
#' ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#' assay(ft2[[2]])
#' aggcounts(ft2[[2]])
#'
#' ## The rowData .n variable tallies number of initial rows that
#' ## were aggregated (irrespective of NAs) for all the samples.
#' rowData(ft2[[2]])
#'
#' ## Ignored when setting na.rm = TRUE
#' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums, na.rm = TRUE)
#' assay(ft3[[2]])
#' aggcounts(ft3[[2]])
#'
#' ## -----------------------------------------------
#' ## Aggregation with missing values in the row data
#' ## -----------------------------------------------
#' ## Row data results without any NAs, which includes the
#' ## Y variables
#' rowData(ft2[[2]])
#'
#' ## Missing value in the Y feature variable
#' rowData(ft_na[[1]])[1, "Y"] <- NA
#' rowData(ft_na[[1]])
#'
#' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
#' ## The Y feature variable has been dropped!
#' assay(ft3[[2]])
#' rowData(ft3[[2]])
#'
#' ## --------------------------------------------
#' ## Using a peptide-by-proteins adjacency matrix
#' ## --------------------------------------------
#'
#' ## Let's use assay peptides from object feat1 and
#' ## define that peptide SYGFNAAR maps to proteins
#' ## Prot A and B
#'
#' se <- feat1[["peptides"]]
#' rowData(se)$Protein[3] <- c("ProtA;ProtB")
#' rowData(se)
#'
#' ## This can also be defined using anadjacency matrix, manual
#' ## encoding here. See PSMatch::makeAdjacencyMatrix() for a
#' ## function that does it automatically.
#' adj <- matrix(0, nrow = 3, ncol = 2,
#'               dimnames = list(rownames(se),
#'                               c("ProtA", "ProtB")))
#' adj[1, 1] <- adj[2, 2] <- adj[3, 1:2] <- 1
#' adj
#'
#' adjacencyMatrix(se) <- adj
#' rowData(se)
#' adjacencyMatrix(se)
#'
#' ## Aggregation using the adjacency matrix
#' se2 <- aggregateFeatures(se, fcol = "adjacencyMatrix",
#'                          fun = MsCoreUtils::colMeansMat)
#'
#' ## Peptide SYGFNAAR was taken into account in both ProtA and ProtB
#' ## aggregations.
#' assay(se2)
#'
#'
#' ## Aggregation by matrix on a QFeature object works as with a
#' ## vector
#' ft <- QFeatures(list(peps = se))
#' ft <- aggregateFeatures(ft, "peps", "adjacencyMatrix", name = "protsByMat",
#'                         fun = MsCoreUtils::colMeansMat)
#' assay(ft[[2]])
#' rowData(ft[[2]])
NULL

#' @exportMethod aggregateQmetadata
#' @rdname qMetadata-aggregate
#' @export
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




#' @export
#'
#' @importFrom ProtGenerics adjacencyMatrix
#'
#' @rdname qMetadata-aggregate
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param qMetaName `character(1)` with the variable name containing
#'     the adjacency matrix. Default is `"qMetadata"`.
#'
#' @param i The index or name of the assays to extract the advaceny
#'     matrix from. All must have a rowdata variable named `qMetaName`.
setMethod("qMetadata", "QFeatures",
          function(object, i, qMetaName = "qMetadata")
            List(lapply(experiments(object)[i],
                        .qMetadata,
                        qMetaName = qMetaName)))

setMethod("qMetadata", "SummarizedExperiment",
          function(object, qMetaName = "qMetadata")
            .qMetadata(object, qMetaName))

#' @export
#'
#' @rdname qMetadata-aggregate
#'
#' @param i When adding an adjacency matrix to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     adjacency matrix will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value An adjacency matrix with row and column names. The
#'     matrix will be coerced to compressed, column-oriented sparse
#'     matrix (class `dgCMatrix`) as defined in the `Matrix` package,
#'     as generaled by the [sparseMatrix()] constructor.
"qMetadata<-" <- function(object, i, qMetaName = "qMetadata", value) {
  if (is.null(colnames(value)) | is.null(rownames(value)))
    stop("The matrix must have row and column names.")
  ## Coerse to a data.frame
  value <- as(value, "data.frame")
  if (inherits(object, "SummarizedExperiment")) {
    if (!identical(rownames(value), rownames(object)))
      stop("Row names of the SummarizedExperiment and the qMetadata data.frame must match.")
    if (qMetaName %in% colnames(rowData(object)))
      stop("Found an existing variable ", qMetaName, ".")
    rowData(object)[[qMetaName]] <- value
    return(object)
  }
  stopifnot(inherits(object, "QFeatures"))
  if (length(i) != 1)
    stop("'i' must be of length one. Repeat the call to add a matrix to multiple assays.")
  if (is.numeric(i) && i > length(object))
    stop("Subscript is out of bounds.")
  if (is.character(i) && !(i %in% names(object)))
    stop("Assay '", i, "' not found.")
  se <- object[[i]]
  object[[i]] <- qMetadata(se, qMetaName) <- value
  return(object)
}

.qMetadata <- function(x, qMetaName = "qMetadata") {
  stopifnot(qMetaName %in% names(rowData(x)))
  ans <- rowData(x)[[qMetaName]]
  if (is.null(colnames(ans)) | is.null(rownames(ans)))
    warning("The qMetadata data.frame should have row and column names.")
  ans
}

# makeAdjacencyMatrix_4Prostar <- function(object, i, fcol, replace = TRUE,...){
#   
#   
#   X <- makeAdjacencyMatrix(rowData(feat2[[2]])[,'Protein_group_IDs'])
#   rownames(X) <- rownames(feat2[[2]])
#   
#   X <- CustomAdjMat(X, mode = 'onlySpec')
#   X <- CustomAdjMat(X, mode = 'topn', n = 2, qData = assay(feat2[[2]]))
#   rowData(feat2[[2]])[['adjacencyMatrix']] <- NULL
#   adjacencyMatrix(feat2[[2]]) <- X
#   
#   
#   X <- makeAdjacencyMatrix(rowData(object[[i]])[,fcol])
#   rownames(X) <- rownames(object[[i]])
#   if (mode == 'topn')
#     X <- CustomAdjMat(X, mode = mode, qData = assay(object[[i]]), n = n)
#   
#   rowData(object[[i]])[['adjacencyMatrix']] <- NULL
#   adjacencyMatrix(object[[i]]) <- CustomAdjMat(object[[i]],
#                                                X, 
#                                                ...)
#   object
#   
# }


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
#' @param ... xxx
#' 
#' @export
#' 
#' @examples
#' feat2 <- readRDS('~/GitHub/DaparToolshedData/data/Exp2_R100_pept.rda')
#' feat2 <- feat2[1:10,]
#' X <- makeAdjacencyMatrix(rowData(feat2[['original_log']])[,'Protein_group_IDs'])
#' rownames(X) <- rownames(feat2[['original_log']])
#' updateAdjacencyMatrix(feat2[['original_log']], X, mode = 'all')
#'
updateAdjacencyMatrix <- function(X, mode = 'all', ...){
  
  X.binary <- X
  X.binary[which(X.binary != 0)] <- 1
  
  switch(mode,
         all = X,
         onlySpec = X[which(rowSums(X.binary) > 1),] <- 0,
         onlyShared = X[which(rowSums(X.binary) == 1),] <- 0,
         topn = .BuildTopnMat(X, ...)
         )
  return(X)
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
.BuildTopnMat <- function(X, qData = NULL, fun = 'rowMedians', n = 10){
  
  xmed <- NULL
  if (fun %in% c('rowMedians', 'rowMeans', 'rowSums')){
    xmed <- as(X * do.call(fun, list(qData)), "dgCMatrix")
  }
  
  # Get the 'n' entities with the best score for each column
  for (c in seq_len(ncol(X))){
    v <- order(xmed[,c],decreasing=TRUE)[seq_len(n)]
    l <- v[which((xmed[,c])[v] != 0)]
    
    if (length(l) > 0){
      diff <- setdiff( which(X[,c] == 1), l)
      if (length(diff)) {X[diff,c] <- 0}
    }
  }
  
  X
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
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
#' obj <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' X <- rowData(obj[[2]])$adjacencyMatrix
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
  #browser()
  X <- as(X, 'matrix')
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' n <- inner.sum(assay(obj[[2]]), X)
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' inner.mean(assay(obj[[2]]), X)
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' inner.aggregate.topn(assay(obj[[2]]), X, n=3)
#' 
#' @importFrom stats median
#' 
inner.aggregate.topn <-function(qData, X, method = 'Mean', n = 10){
  
  X <- as.matrix(X)
  
  med <- apply(qData, 1, median)
  xmed <- as(X * med, "dgCMatrix")
  for (c in seq_len(ncol(X))){
    v <- order(xmed[,c],decreasing=TRUE)[seq_len(n)]
    l <- v[which((xmed[,c])[v] != 0)]
    
    if (length(l) > 0){
      diff <- setdiff( which(X[,c] == 1), l)
      if (length(diff)) {X[diff,c] <- 0}
    }
  }
  
  Mp <- NULL
  switch(method,
         Mean = Mp <- inner.mean(qData, X),
         Sum = Mp <- inner.sum(qData, X)
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- GetAdjMat(obj[[2]])$all
#' qPepData <- assay(obj[[2]])
#' inner.aggregate.iter(qPepData, X)
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
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' obj <- addListAdjacencyMatrices(obj, 2)
#' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' aggTopn(assay(obj,2), X, n=3)
#' 
#' @export
#' 
aggTopn <- function(qData, X,  method = 'Mean', n = 10){
  inner.aggregate.topn(qData, X, method, n)
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



#-----------------------------------------

#' @title
#' Combine peptide metadata to build protein metadata
#' 
#' @description 
#' Agregation rules for the cells quantitative metadata of peptides. 
#' Please refer to the qMetadata.def vocabulary in `qMetadata.def()`
#' 
#' # Basic agreagtion
#' Agregation of non imputed values (2.X) with quantitative values 
#' (1.0, 1.X, 3.0, 3.X)
#' |----------------------------
#' Not possible
#' |----------------------------
#' 
#' Agregation of different types of missing values (among 2.1, 2.2)
#' |----------------------------
#' * Agregation of 2.1 peptides between each other gives a missing value 
#'   non imputed (2.0)
#' * Agreagtion of 2.2 peptides between each other givesa missing value 
#'   non imputed (2.0)
#' * Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed (2.0)
#' |----------------------------
#' 
#' 
#' Agregation of a mix of quantitative values (among 1.0, 1.1, 1.2, 3.0, 3.X)
#' |----------------------------
#' * if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, 
#'   then the final metadata is set the this tag
#' * if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, 
#'   then the final metadata is set to 1.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0, 
#'   then the final metadata is set to 3.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0 and other (X.X),
#'   then the final metadata is set to 4.0
#' |----------------------------
#' 
#' # Post processing
#' Update metacell with POV/MEC status for the categories 2.0 and 3.0
#' TODO
#' 
#' @param met xxx
#' 
#' @param level xxx
#' 
#' @examples
#' \dontrun{
#' ll <- qMetadata.def('peptide')$node
#' for (i in 1:length(ll))
#' test <- lapply(combn(ll, i, simplify = FALSE), 
#' function(x) tag <- qMetadata_combine(x, 'peptide'))
#' }
#' 
#' 
qMetadata_combine <- function(met, level) {
  tag <- NULL
  if (length(met)==0)
    return('missing')
  
  u_met <- unique(met)
  
  # Define an auxiliary function
  ComputeNbTags <- function(tag){
    sum(unlist(lapply( search.qMetadata.tags(tag, level), 
                       function(x) length(grep(x, u_met)))))
  }
  
  
  nb.tags <- lapply(qMetadata.def(level)$node, 
                    function(x) as.numeric(x %in% u_met))
  n.imputed <- ComputeNbTags('imputed')
  n.missing <- ComputeNbTags('missing')
  n.quanti <- ComputeNbTags('quanti')
  
  
  if(n.missing > 0 && (n.imputed > 0 || n.quanti > 0)) tag <- 'STOP'
  # stop("You try to combine missing values (2.X) with quantitative values (1.X or 3.X).")
  
  # sw : Agregation of a mix of 2.X gives a missing value non imputed (2.0)
  if (n.missing > 0 && n.quanti == 0 && n.imputed == 0) tag <- 'missing'
  
  
  # # Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed (2.0)
  # if (length(u_met)== length(grep('missing_', u_met))) tag <- 'missing'
  # 
  # # Agreagtion of 2.2 peptides between each other givesa missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_MEC' == u_met) tag <- 'missing'
  # 
  # # Agreagtion of 2.2 peptides between each other gives a missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_POV' == u_met) tag <- 'missing'
  #     
  # # Agregation of 2.1 peptides between each other gives a missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_MEC' == u_met) tag <- 'missing'
  
  # if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, then the final
  # metadata is set the this tag
  if (length(u_met)==1 && u_met == 'quanti') tag <- 'quanti'
  if (length(u_met)==1 && u_met == 'identified') tag <- 'identified'
  if (length(u_met)==1 && u_met == 'recovered') tag <- 'recovered'
  
  
  # if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, then the final
  # metadata is set to 1.0
  if (n.quanti > 1 && n.imputed == 0 && n.missing==0) tag <- 'quanti'
  
  
  # If the set of metacell to agregate is a mix of 3.X and 3.0, then the final
  # metadata is set to 3.0
  if (n.quanti == 0 && n.imputed > 0 && n.missing == 0) tag <- 'imputed'
  
  # If the set of metacell to agregate is a mix of 3.X and 3.0 and other (X.X), 
  # then the final metadata is set to 4.0
  if (n.quanti > 0 && n.imputed > 0 && n.missing == 0)
    tag <- 'combined'
  
  #print(paste0(paste0(u_met, collapse=' '), ' ---> ', tag))
  return(tag)
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
#' @rdname qMetadata-aggregate
#' 
aggQmeta <- function(object, conds) {
  stopifnot(inherits(object, "SummarizedExperiment"))
  
  qMeta = qMetadata(object)
  level = GetTypeDataset(object)
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
#' @rdname qMetadata-aggregate
#' 
FinalizeAggregation <- function(object){
  
  object.from <- object[[from]]
  object.to <- object[[to]]
  
  fromData <- assay(object.from)
  toData <- assay(object.to)
  
  
  #(obj.pep, pepData, protData, metacell, X)
  
  protData <- as.matrix(protData)
  X <- as.matrix(X)
  protData[protData==0] <- NA
  protData[is.nan(protData)] <- NA
  protData[is.infinite(protData)] <-NA
  
  temp <- GetDetailedNbPeptidesUsed(X, pepData)
  
  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", colnames(pepData), sep="")
  rownames(pepSharedUsed) <- colnames(X)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", colnames(pepData), sep="")
  rownames(pepSpecUsed) <- colnames(X)
  
  pepTotalUsed <- as.matrix(GetNbPeptidesUsed(X, pepData))
  colnames(pepTotalUsed) <- paste("pepTotal.used.", colnames(pepData), sep="")
  rownames(pepTotalUsed) <- colnames(X)
  
  n <- GetDetailedNbPeptides(X)
  
  
  fd <- data.frame(proteinId = rownames(protData),
                   nPepTotal = n$nTotal,
                   nPepShared = n$nShared, 
                   nPepSpec = n$nSpec, 
                   pepSpecUsed, 
                   pepSharedUsed, 
                   pepTotalUsed, 
                   protMetacell)
  rownames(fd) <- colnames(X)
  obj.prot <- MSnSet(exprs = log2(protData), 
                     fData = fd, 
                     pData = Biobase::pData(obj.pep))
  
  obj.prot@experimentData@other <- obj.pep@experimentData@other
  obj.prot@experimentData@other$typeOfData <- "protein"
  obj.prot@experimentData@other$keyId <- 'proteinId'
  obj.prot@experimentData@other$proteinId <- 'proteinId'
  
  
  obj.prot@experimentData@other$Prostar_Version <- NA
  obj.prot@experimentData@other$DAPAR_Version <- NA
  
  tryCatch({
    find.package("Prostar")
    find.package("DaparToolshed")
    
    obj.prot@experimentData@other$Prostar_Version <- Biobase::package.version('Prostar')
    obj.prot@experimentData@other$DAPAR_Version <- Biobase::package.version('DaparToolshed')
  },
  error = function(e) {
    obj.prot@experimentData@other$Prostar_Version <- NA
    obj.prot@experimentData@other$DaparToolshed_Version <- NA
  }
  )
  
  return (obj.prot)
}




