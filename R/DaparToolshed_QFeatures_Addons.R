##' @title 
##' 
##' Additional accessors for instances of class `Qfeatures`.
##' 
##' @description
##'  
##' These names are common to all assays contained in the object. This is why
##' they are stored in the global metadata. This function is used whenever it is necessary
##' to (re)detect MEC and POV (new dataset or when post processing protein qMetadata 
##' after aggregation)
##' 
##' @param object An instance of class `SummarizedExperiment` or `QFeatures`.
##' 
##' @param i The index or name of the assays to extract the quantitative metadata
##'  from. All must have a rowdata variable named as `slotName`
##' 
##' @param slotName xxx
##' 
##' @param value xxx
##' 
##' 
##' @details 
##' 
##' Additional slots for rowdata of a `SummarizedExperiment` object:
##'  - qMetadata: xxx
##' 
##' Additional slots for Metadata for a `QFeatures` object:
##'  - xxx: xxxx
##' 
##' Additional slots for Metadata for a `SummarizedExperiment` object:
##'  - qMetadata: xxxx
##'  - parentProtId: xxx
##'  - idcol: xxxx
##'  - typeDataset: xxx
##'  
##' 
##' @section Quantitative metadata:
##' 
##' Default slotName is `"qMetadata"`.
##'  The value is an adjacency matrix with row and column names. The
##'     matrix will be coerced to compressed, column-oriented sparse
##'     matrix (class `dgCMatrix`) as defined in the `Matrix` package,
##'     as generaled by the [sparseMatrix()] constructor.
##' 
##' @name QFeatures-accessors
##' @rdname QFeatures-accessors
##'   
##' @examples 
##' data(ft)
##' design(ft)
NULL

##' @exportMethod qMetadata
##' @rdname QFeatures-accessors
setMethod("qMetadata", "QFeatures",
          function(object, i, slotName = "qMetadata")
            List(lapply(experiments(object)[i],
                        function(x)
                          .GetRowdataSlot(x, slotName = slotName)))
)




##' @export
##' @rdname QFeatures-accessors
setMethod("qMetadata", "SummarizedExperiment",
          function(object, slotName = "qMetadata")
            .GetRowdataSlot(object, slotName))



##' @export
##' @rdname QFeatures-accessors
"qMetadata<-" <- function(object,i,slotName = "qMetadata", value) {
  if (is.null(colnames(value)) | is.null(rownames(value)))
    stop("The DataFrame must have row and column names.")
  ## Coerse to a data.frame
  #value <- as(value, "data.frame")
  if (inherits(object, "SummarizedExperiment")) {
    if (!identical(rownames(value), rownames(object)))
      stop("Row names of the SummarizedExperiment and the DataFrame must match.")
    #if (slotName %in% colnames(rowData(object)))
    #  stop("Found an existing variable ", slotName, ".")
    rowData(object)[[slotName]] <- value
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
  object[[i]] <- qMetadata(se, slotName) <- value
  return(object)
}

##' @title .GetMetadataSlot
##' @noRd
.GetMetadataSlot <- function(x, slotName = NULL) {
  ans <- metadata(x)[[slotName]]
  ans
}

##' @title .GetMetadataSlot
##' @noRd
.GetRowdataSlot <- function(x, slotName = NULL) {
  ans <- rowData(x)[[slotName]]
  ans
}


##' @exportMethod typeDataset
##' @rdname QFeatures-accessors
setMethod("typeDataset", "QFeatures",
          function(object, i, slotName = "typeDataset")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)))
##' @export
##' @rdname QFeatures-accessors
setMethod("typeDataset", "SummarizedExperiment",
          function(object, slotName = "typeDataset")
            .GetMetadataSlot(object, slotName))


##' @export
##' @rdname QFeatures-accessors
"typeDataset<-" <- function(object, i, slotName = "typeDataset", value) {
  if (inherits(object, "SummarizedExperiment")) {
    metadata(object)[[slotName]] <- value
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
  metadata(object[[i]])[[slotName]] <- value
  return(object)
}



##' @exportMethod idcol
##' @rdname QFeatures-accessors
setMethod("idcol", "QFeatures",
          function(object, i, slotName = "idcol")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)))
##' @export
##' @rdname QFeatures-accessors
setMethod("idcol", "SummarizedExperiment",
          function(object, slotName = "idcol")
            .GetMetadataSlot(object, slotName))


##' @export
##' @rdname QFeatures-accessors
"idcol<-" <- function(object, i, slotName = "idcol", value) {
  if (inherits(object, "SummarizedExperiment")) {
    metadata(object)[[slotName]] <- value
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
  metadata(object[[i]])[[slotName]] <- value
  return(object)
}




##' @exportMethod parentProtId
##' @rdname QFeatures-accessors
setMethod("parentProtId", "QFeatures",
          function(object, i, slotName = "parentProtId")
            List(lapply(experiments(object)[i],
                        .parentProtId,
                        slotName = slotName)))

setMethod("parentProtId", "SummarizedExperiment",
          function(object, slotName = "parentProtId")
            if (typeDataset(object) == 'peptide')
              .GetMetadataSlot(object, slotName))


##' @export
##' @rdname QFeatures-accessors
"parentProtId<-" <- function(object, i, slotName = "parentProtId", value) {
  if (inherits(object, "SummarizedExperiment")) {
    if (typeDataset(object) != 'peptide')
      stop("The dataset must contain peptides.")
    metadata(object)[[slotName]] <- value
    return(object)
  }
  stopifnot(inherits(object, "QFeatures"))
  if (typeDataset(object[[i]]) != 'peptide')
    stop("The dataset must contain peptides.")
  if (length(i) != 1)
    stop("'i' must be of length one. Repeat the call to add a matrix to multiple assays.")
  if (is.numeric(i) && i > length(object))
    stop("Subscript is out of bounds.")
  if (is.character(i) && !(i %in% names(object)))
    stop("Assay '", i, "' not found.")
  se <- object[[i]]
  metadata(object[[i]])[[slotName]] <- value
  return(object)
}





##' @exportMethod analysis
##' @rdname QFeatures-accessors
setMethod("analysis", "QFeatures",
          function(object, i, slotName = "analysis")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)
            ))
##' @export
##' @rdname QFeatures-accessors
setMethod("analysis", "SummarizedExperiment",
          function(object, slotName = "analysis")
            .GetMetadataSlot(object, slotName))


##' @export
##' @rdname QFeatures-accessors
"analysis<-" <- function(object, i, slotName = "analysis", value) {
  if (inherits(object, "SummarizedExperiment")) {
    metadata(object)[[slotName]] <- value
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
  metadata(object[[i]])[[slotName]] <- value
  return(object)
}





##' @exportMethod version
##' @rdname QFeatures-accessors
setMethod("version", "QFeatures",
          function(object, slotName = "version")
            .GetMetadataSlot(object, slotName = slotName)
)



##' @export
##' @rdname QFeatures-accessors
"version<-" <- function(object, slotName = "version", value) {
  stopifnot (inherits(object, "QFeatures"))
  metadata(object)[[slotName]] <- value
  return(object)
}


##' @exportMethod design
##' @rdname QFeatures-accessors
setMethod("design", "QFeatures",
          function(object, slotName = "design")
            colData(object)
)



##' @export
##' @rdname QFeatures-accessors
"design<-" <- function(object, slotName = "design", value) {
  stopifnot (inherits(object, "QFeatures"))
  colData(object)@listData <- value
  return(object)
}

