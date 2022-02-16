
#' @title Additional accessors to the quantitative cell metadata
#' 
#' @description 
#' These names are common to all assays contained in the object. This is why
#' they are stored in the global metadata. This function is used whenever it is necessary
#' to (re)detect MEC and POV (new dataset or when post processing protein qMetadata 
#' after aggregation)
#' 
#' @details 
#' 
#' 
#' Additional slots for rowdata of a `SummarizedExperiment` object
#'  - qMetadata: xxx
#' 
#' Additional slots for Metadata for a `QFeatures` object
#'  - xxx: xxxx
#' 
#' Additional slots for Metadata for a `SummarizedExperiment` object
#'  - qMetadata: xxxx
#'  - parentProtId: xxx
#'  - idcol: xxxx
#'  - typeDataset: xxx
#'  
#'  
NULL

#--------------------------------

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the quantitative metadata. Default is `"qMetadata"`.
#'
#' @param i The index or name of the assays to extract the quantitative metadata
#'  from. All must have a rowdata variable named `qMetadata`.
setMethod("qMetadata", "QFeatures",
          function(object, i, slotName = "qMetadata")
            List(lapply(experiments(object)[i],
                        .qMetadata,
                        slotName = slotName)))

setMethod("qMetadata", "SummarizedExperiment",
          function(object, slotName = "qMetadata")
            .qMetadata(object, slotName))


#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an adjacency matrix to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     adjacency matrix will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value An `DataFrame` with row and column names.
"qMetadata<-" <- function(object, i, slotName = "qMetadata", value) {
  if (is.null(colnames(value)) | is.null(rownames(value)))
    stop("The DataFrame must have row and column names.")
  
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

.qMetadata <- function(x, slotName = "qMetadata") {
  #stopifnot(slotName %in% names(rowData(x)))
  ans <- rowData(x)[[slotName]]
  if (is.null(colnames(ans)) | is.null(rownames(ans)))
    warning("The qMetadata dataframe should have row and column names.")
  ans
}


#--------------------------------------------------------


#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the adjacency matrix. Default is `"typeDataset"`.
#'
#' @param i The index or name of the assays to extract the advaceny
#'     matrix from. All must have a rowdata variable named `slotName`.
setMethod("typeDataset", "QFeatures",
          function(object, i, slotName = "typeDataset")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)))

setMethod("typeDataset", "SummarizedExperiment",
          function(object, slotName = "typeDataset")
            .GetMetadataSlot(object, slotName))

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an type dataset to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     type dataset  will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value A `character()`.
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


#------------------------------------------------------------



#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the key Id Column of the se. Default is `"idcol"`.
#'
#' @param i The index or name of the assays to extract the information from.
#'  All must have a rowdata variable named `slotName`.
setMethod("idcol", "QFeatures",
          function(object, i, slotName = "idcol")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)))

setMethod("idcol", "SummarizedExperiment",
          function(object, slotName = "idcol")
            .GetMetadataSlot(object, slotName))

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an type dataset to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     type dataset  will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value A `character()`.
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



#--------------------------


#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the name of the column containing. Default is `"parentProtId"`.
#'
#' @param i The index or name of the assays to extract the parentProtId from. 
#'          All must have a rowdata variable named `slotName`.
setMethod("parentProtId", "QFeatures",
          function(object, i, slotName = "parentProtId")
            List(lapply(experiments(object)[i],
                        .parentProtId,
                        slotName = slotName)))

setMethod("parentProtId", "SummarizedExperiment",
          function(object, slotName = "parentProtId")
            if (typeDataset(object) == 'peptide')
              .GetMetadataSlot(object, slotName))

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an type dataset to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     type dataset  will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value A `character(1)`.
#' 
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




#--------------------------


#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the name of the column containing. Default is `"analysis"`.
#'
#' @param i The index or name of the assays to extract the analysis from. 
#'          All must have a metadata variable named `slotName`.
setMethod("analysis", "QFeatures",
          function(object, i, slotName = "analysis")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)
                        ))

setMethod("analysis", "SummarizedExperiment",
          function(object, slotName = "analysis")
            .GetMetadataSlot(object, slotName))

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an type dataset to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     type dataset  will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value A `character(1)`.
#' 
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



#--------------------------


#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or
#'     `QFeatures`.
#'
#' @param slotName `character(1)` with the variable name containing
#'     the versions of packages used. Default is `"version"`.
#'
#' @param i The index or name of the assays to extract the version from. 
#'          All must have a rowdata variable named `slotName`.
setMethod("version", "QFeatures",
          function(object, i, slotName = "version")
            List(lapply(experiments(object)[i],
                        .GetMetadataSlot,
                        slotName = slotName)))

setMethod("version", "SummarizedExperiment",
          function(object, slotName = "version")
            .GetMetadataSlot(object, slotName))

#' @export
#'
#' @rdname QFeatures-supplementary-accessors
#'
#' @param i When adding an type dataset to an assay of a
#'     `QFeatures` object, the index or name of the assay the
#'     type dataset  will be added to. Ignored when `x` is an
#'     `SummarizedExperiment`.
#'
#' @param value A `character(1)`.
#' 
"version<-" <- function(object, i, slotName = "version", value) {
  if (inherits(object, "SummarizedExperiment")) {
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


.GetMetadataSlot <- function(x, slotName = "version") {
  ans <- metadata(x)[[slotName]]
  ans
}