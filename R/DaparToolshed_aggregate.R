

#' @title Aggregate an assay of peptides into proteins
#' 
#' @description 
#' 
#' This function is a wrapper for the function [QFeatures::aggregateFeatures()].
#' It just add some meta informations to the aggregated assay, specific to
#' DaparToolshed
#' 
##' @title Aggregate an assay's quantitative features
##'
##' @description
##'
##' This function aggregates the quantitative features of an assay,
##' applying a summarisation function (`fun`) to sets of features.
##' The `fcol` variable name points to a rowData column that defines
##' how to group the features during aggregate. This variable can
##' eigher be a vector (we then refer to an *aggregation by vector*)
##' or an adjacency matrix (*aggregation by matrix*).
##'
##' The rowData of the aggregated `SummarizedExperiment` assay
##' contains a `.n` variable that provides the number of parent
##' features that were aggregated.
##'
##' When aggregating with a vector, the newly aggregated
##' `SummarizedExperiment` assay also contains a new `aggcounts` assay
##' containing the aggregation counts matrix, i.e. the number of
##' features that were aggregated for each sample, which can be
##' accessed with the `aggcounts()` accessor.
##'
##' @param object An instance of class [QFeatures] or [SummarizedExperiment].
##'
##' @param i The index or name of the assay which features will be
##'     aggregated the create the new assay.
##'
##' @param fcol A `character(1)` naming a rowdata variable (of assay
##'     `i` in case of a `QFeatures`) defining how to aggregate the
##'     features of the assay. This variable is either a `character`
##'     or a (possibly sparse) matrix. See below for details.
##'
##' @param name A `character(1)` naming the new assay. Default is
##'     `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##'
##' @param fun A function used for quantitative feature
##'     aggregation. See Details for examples.
##'
##' @param ... Additional parameters passed the `fun`.
##'
##' @return A `QFeatures` object with an additional assay or a
##'  `SummarizedExperiment` object (or subclass thereof).
##'
##' @details
##'
##' Aggregation is performed by a function that takes a matrix as
##' input and returns a vector of length equal to `ncol(x)`. Examples
##' thereof are
##'
##' - [MsCoreUtils::medianPolish()] to fits an additive model (two way
##'   decomposition) using Tukey's median polish_ procedure using
##'   [stats::medpolish()];
##'
##' - [MsCoreUtils::robustSummary()] to calculate a robust aggregation
##'   using [MASS::rlm()] (default);
##'
##' - [base::colMeans()] to use the mean of each column;
##'
##' - `colMeansMat(x, MAT)` to aggregate feature by the calculating
##'    the mean of peptide intensities via an adjacency matrix. Shared
##'    peptides are re-used multiple times.
##'
##' - [matrixStats::colMedians()] to use the median of each column.
##'
##' - [base::colSums()] to use the sum of each column;
##'
##' - `colSumsMat(x, MAT)` to aggregate feature by the summing the
##'    peptide intensities for each protein via an adjacency
##'    matrix. Shared peptides are re-used multiple times.
##'
##' See [MsCoreUtils::aggregate_by_vector()] for more aggregation functions.
##'
##' @section Missing quantitative values:
##'
##' Missing quantitative values have different effects based on the
##' aggregation method employed:
##'
##' - The aggregation functions should be able to deal with missing
##'   values by either ignoring or propagating them. This is often
##'   done with an `na.rm` argument, that can be passed with
##'   `...`. For example, `rowSums`, `rowMeans`, `rowMedians`,
##'   ... will ignore `NA` values with `na.rm = TRUE`, as illustrated
##'   below.
##'
##' - Missing values will result in an error when using `medpolish`,
##'   unless `na.rm = TRUE` is used. Note that this option relies on
##'   implicit assumptions and/or performes an implicit imputation:
##'   when summing, the values are implicitly imputed by 0, assuming
##'   that the `NA` represent a trully absent features; when
##'   averaging, the assumption is that the `NA` represented a
##'   genuinely missing value.
##'
##' - When using robust summarisation, individual missing values are
##'   excluded prior to fitting the linear model by robust
##'   regression. To remove all values in the feature containing the
##'   missing values, use [filterNA()].
##'
##' More generally, missing values often need dedicated handling such
##' as filtering (see [filterNA()]) or imputation (see [impute()]).
##'
##' @section Missing values in the row data:
##'
##' Missing values in the row data of an assay will also impact the
##' resulting (aggregated) assay row data, as illustrated in the
##' example below. Any feature variables (a column in the row data)
##' containing `NA` values will be dropped from the aggregated row
##' data. The reasons underlying this drop are detailed in the
##' `reduceDataFrame()` manual page: only invariant aggregated rows,
##' i.e. rows resulting from the aggregation from identical variables,
##' are preserved during aggregations.
##'
##' The situation illustrated below should however only happen in rare
##' cases and should often be imputable using the value of the other
##' aggregation rows before aggregation to preserve the invariant
##' nature of that column. In cases where an `NA` is present in an
##' otherwise variant column, the column would be dropped anyway.
##'
##' @section Using an adjacency matrix:
##'
##' When considering non-unique peptides explicitly, i.e. peptides
##' that map to multiple proteins rather than as a protein group, it
##' is convenient to encode this ambiguity explicitly using a
##' peptide-by-proteins (sparse) adjacency matrix. This matrix is
##' typically stored in the rowdata and set/retrieved with the
##' [adjacencyMatrix()] function. It can be created manually (as
##' illustrated below) or using `PSMatch::makeAdjacencyMatrix()`.
##'
##' @seealso The *QFeatures* vignette provides an extended example and
##'     the *Processing* vignette, for a complete quantitative
##'     proteomics data processing pipeline. The
##'     [MsCoreUtils::aggregate_by_vector()] manual page provides
##'     further details.
##'
##' @aliases aggregateFeatures aggregateFeatures,QFeatures-method
##'     aggcounts aggcounts,SummarizedExperiment-method
##'     adjacencyMatrix,SummarizedExperiment-method
##'     adjacencyMatrix,QFeatures-method
##'
##' @name aggregateFeatures
##'
##' @rdname QFeatures-aggregate
##'
##' @importFrom MsCoreUtils aggregate_by_vector aggregate_by_matrix robustSummary colCounts
##'
##' @examples
##'
##' ## ---------------------------------------
##' ## An example QFeatures with PSM-level data
##' ## ---------------------------------------
##' data(feat1)
##' feat1
##'
##' ## Aggregate PSMs into peptides
##' feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
##' feat1
##'
##' ## Aggregate peptides into proteins
##' feat1 <- aggregateFeatures(feat1, "peptides", "Protein", name = "proteins")
##' feat1
##'
##' assay(feat1[[1]])
##' assay(feat1[[2]])
##' aggcounts(feat1[[2]])
##' assay(feat1[[3]])
##' aggcounts(feat1[[3]])
##'
##' ## --------------------------------------------
##' ## Aggregation with missing quantitative values
##' ## --------------------------------------------
##' data(ft_na)
##' ft_na
##'
##' assay(ft_na[[1]])
##' rowData(ft_na[[1]])
##'
##' ## By default, missing values are propagated
##' ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
##' assay(ft2[[2]])
##' aggcounts(ft2[[2]])
##'
##' ## The rowData .n variable tallies number of initial rows that
##' ## were aggregated (irrespective of NAs) for all the samples.
##' rowData(ft2[[2]])
##'
##' ## Ignored when setting na.rm = TRUE
##' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums, na.rm = TRUE)
##' assay(ft3[[2]])
##' aggcounts(ft3[[2]])
##'
##' ## -----------------------------------------------
##' ## Aggregation with missing values in the row data
##' ## -----------------------------------------------
##' ## Row data results without any NAs, which includes the
##' ## Y variables
##' rowData(ft2[[2]])
##'
##' ## Missing value in the Y feature variable
##' rowData(ft_na[[1]])[1, "Y"] <- NA
##' rowData(ft_na[[1]])
##'
##' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
##' ## The Y feature variable has been dropped!
##' assay(ft3[[2]])
##' rowData(ft3[[2]])
##'
##' ## --------------------------------------------
##' ## Using a peptide-by-proteins adjacency matrix
##' ## --------------------------------------------
##'
##' ## Let's use assay peptides from object feat1 and
##' ## define that peptide SYGFNAAR maps to proteins
##' ## Prot A and B
##'
##' se <- feat1[["peptides"]]
##' rowData(se)$Protein[3] <- c("ProtA;ProtB")
##' rowData(se)
##'
##' ## This can also be defined using anadjacency matrix, manual
##' ## encoding here. See PSMatch::makeAdjacencyMatrix() for a
##' ## function that does it automatically.
##' adj <- matrix(0, nrow = 3, ncol = 2,
##'               dimnames = list(rownames(se),
##'                               c("ProtA", "ProtB")))
##' adj[1, 1] <- adj[2, 2] <- adj[3, 1:2] <- 1
##' adj
##'
##' adjacencyMatrix(se) <- adj
##' rowData(se)
##' adjacencyMatrix(se)
##'
##' ## Aggregation using the adjacency matrix
##' se2 <- aggregateFeatures(se, fcol = "adjacencyMatrix",
##'                          fun = MsCoreUtils::colMeansMat)
##'
##' ## Peptide SYGFNAAR was taken into account in both ProtA and ProtB
##' ## aggregations.
##' assay(se2)
##'
##'
##' ## Aggregation by matrix on a QFeature object works as with a
##' ## vector
##' ft <- QFeatures(list(peps = se))
##' ft <- aggregateFeatures(ft, "peps", "adjacencyMatrix", name = "protsByMat",
##'                         fun = MsCoreUtils::colMeansMat)
##' assay(ft[[2]])
##' rowData(ft[[2]])
NULL

##' @exportMethod aggregateFeatures4Prostar
##' @rdname DaparToolshed-aggregate
setMethod("aggregateFeatures4Prostar", "QFeatures",
          function(object, i, fcol, name = "newAssay",
                   fun = MsCoreUtils::robustSummary, ...) {
            if (isEmpty(object))
              return(object)
            if (name %in% names(object))
              stop("There's already an assay named '", name, "'.")
            if (missing(i))
              i <- main_assay(object)
            
            
            
            # Add stats on agregation
            aggAssay <- aggregateFeatures4Prostar(object[[i]], fcol, fun, ...)
            
            ## Add the assay to the QFeatures object
            object <- addAssay(object,
                               aggAssay,
                               name = name)
            ## Link the input assay to the aggregated assay
            addAssayLink(object,
                         from = i,
                         to  = name,
                         varFrom = fcol,
                         varTo = fcol)
          })


##' @exportMethod aggregateFeatures4Prostar
##' @rdname DaparToolshed-aggregate
setMethod("aggregateFeatures4Prostar", "SummarizedExperiment",
          function(object, fcol, fun = MsCoreUtils::robustSummary, ...)
            .aggregateFeatures4Prostar(object, fcol, fun, ...))


.aggregateFeatures4Prostar <- function(object, fcol, fun, ...) {
  
  ## Create the aggregated assay
  aggAssay <- aggregateFeatures(object, fcol, fun, ...)
  
  
  # add agregation of qMetadata
  # Aggregate the quantitative metdata
  aggQ <- aggQmetadata(object, conds = conds)
  qMetadata(aggAssay) <- aggQ
  
  typeDataset(aggAssay) <- 'proteins'
  idcol(aggAssay) <- NULL

  return(aggAssay)
}



#----------------------------------------



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
#' agg.qMmetadata <- aggQmeta(obj, 2, X)
#' agg.meta <- aggQmeta(obj[[2]], X)
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
aggQmetadata <- function(object, conds) {
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

