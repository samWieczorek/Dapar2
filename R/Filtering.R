
#' @title Complex filter features
#' 
#' @description Filter an assay with a complex query
#' 
#' @author Samuel Wieczorek
#' 
"ComplexFilterFeatures"


#' @param xxx
#' 
#' @export
#' 
setMethod("ComplexFilterFeatures", "SummarizedExperiment",
          function(object, cmd, ...) {
            x <- rowData(object)
            
            # Build fictive column that represents the complex query
            
            # i <- length(experiments(object))
            # 
            # ## Create a fictive column for each assays otherwise the filter method
            # ## of QFeatures will truncate the object
            # for (k in 1:length(experiments(object)))
            #   rowData(object[[k]]) <- setNames(cbind(rowData(object[[k]]),tmp=rep(0, nrow(object[[k]]))),
            #                                    c(names(rowData(object[[k]])),newColName)
            #   )
            if (missing(cmd))
              stop("'cmd' is required;")
            else if (!(cmd %in% c('delete', 'keep')))
              stop("'cmd' must be one of the following values: `delete` or `keep`.")
            
            
            rowData(object) <- cbind(rowData(object),
                                     fictive.col = rep(0, nrow(object))
                                     )
            
            indices <- GetIndices_ComplexQueryFiltering(object, level, pattern, type, percent, op, th )
            
            rowData(object)[indices, 'fictive.col'] <- 1 
            
            # Create VariableFilter
            var.filter <- VariableFilter(field = 'fictive.col', value = cmd=='keep', condition = '==')
            
            # Filter dataset
            object <- filterFeaturesSam(object, i = 2, filter = var.filter)
            
            # Remove fictive column
            rowData(object) <- rowData(object)[, -'fictive.col']
            
            object
          }
)


#' @param  object An object of class `QFeatures`.
#' 
#' @param i A numeric vector or a character vector giving the index or the 
#'     name, respectively, of the assay(s) to be processed.
#'
#' @param name A `character(1)` naming the new assay name. Defaults
#'     are `ttestAssay`.
#' 
#' @param filter xxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @rdname filterFeaturesSam
#' 
#' @export
#' 
setMethod("ComplexFilterFeatures", "QFeatures",
          function(object, i, name = "filterAssay", ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            
            
            
            argg <- c(as.list(environment()))
            tmp <- ComplexFilterFeatures(object[[i]], ...)
            
            if (nrow(tmp) == 0){
              warning('The filtering has not been proceeded beacause it empties all the dataset.')
              ## The object is not affected
              object
            } else {
              object <- addAssay(object,
                                 tmp,
                                 name)
              addAssayLink(object, from = i, to = name)
            }
          }
          
)








#' @title Check the validity of the experimental design
#'
#' @description
#'
#' This manual page describes the The filtering functions of DAPAR have not been all moved to DaparToolshed as we now use the QFeatures package
#' which provides some filtering functions, especially on features that are present in the rowData
#' of the datasets.
#'  The filtering functions on numerical values are deleted because the same functions exist in QFeatures.
#'  For the missing values filtering on conditions, we do not use DAPAR functions anymore. Instead, we use
#'  the numerical filtering functions in QFeatures. To do so, it is necessary to build some rowdata for the
#'  SummarizedExperiment (not necessary stored in the object) which count the number of missing values w.r.t.
#'  the type of filtering: whole matrix, at least one value per condition, etc...
#'
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_prot, package='DAPARdata2')
#' object <- Exp1_R25_prot
#' filter_600 <- VariableFilter(field='Sequence_length', value=as.numeric('600'), condition='>')
#' object <- filterFeaturesSam(object, i=2, filter=reverse_600)
#'
#'
#'
#'
"filterFeaturesSam"

#' @param  object An object of class `QFeatures`.
#'
#' @param filter xxx
#'
#' @param ... Additional parameters passed to inner functions.
#'
#' @export
#'
#' @rdname filterFeaturesSam
#'
#' @importFrom AnnotationFilter field value condition
#' 
setMethod("filterFeaturesSam", "SummarizedExperiment",
          function(object, filter, ...) {
            x <- rowData(object)
            if (field(filter) %in% names(x)){
              sel <- do.call(condition(filter),
                             list(x[, field(filter)], value(filter))
              )
            } else {
              sel <- rep(FALSE, nrow(x))
            }
            object[sel, ]
          }
)


#' @param  object An object of class `QFeatures`.
#' 
#' @param i A numeric vector or a character vector giving the index or the 
#'     name, respectively, of the assay(s) to be processed.
#'
#' @param name A `character(1)` naming the new assay name. Defaults
#'     are `ttestAssay`.
#' 
#' @param filter xxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @rdname filterFeaturesSam
#' 
#' @export
#' 
setMethod("filterFeaturesSam", "QFeatures",
          function(object, i, name = "filterAssay", filter, ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            
            
            
            argg <- c(as.list(environment()))
            tmp <- filterFeaturesSam(object[[i]], filter)
            
            if (nrow(tmp) == 0){
              warning('The filtering has not been proceeded beacause it empties all the dataset.')
              ## The object is not affected
              object
            } else {
              object <- addAssay(object,
                                 tmp,
                                 name)
              addAssayLink(object, from = i, to = name)
            }
          }
          
)


#' @title Quantitative cell metadata scopes for filtering
#' 
#' @export
#' 
MetacellFilteringScope <- function()
  c("None", "WholeLine", "WholeMatrix", "AllCond", "AtLeastOneCond")



#' @title Operators for complex queries
#' 
#' @export
#' 
SymFilteringOperators <- function()
  c('<=','<', '>=', '>', '==', '!=')




#' This function returns the indice of the lines to delete, based on a 
#' prefix string
#' 
#' @title Get the indices of the lines to delete, based on a prefix string
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param idLine2Delete The name of the column that correspond to the data 
#' to filter
#' 
#' @param prefix A character string that is the prefix to find in the data
#' 
#' @return A vector of integers.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' ind <- getIndicesOfLinesToRemove(Exp1_R25_pept[1:100], "Potential_contaminant", 
#' prefix="+")
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
GetIndicesForPrefix <- function(object, col = NULL, prefix = NULL)
{
  if ((prefix == "") || is.null(prefix)) {
    # warning ("No change was made")
    return (NULL)
  }
  
  t <- prefix == substring(rowdata(object)[,col], 1, nchar(prefix))

  return(which( t == TRUE))
}



#' @title Delete the lines in the matrix of intensities and the metadata table
#' given their indice.
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param level A vector of integers which are the indices of lines to 
#' delete.
#' 
#' @param pattern A string to be included in the \code{MSnSet}
#' object for log. 
#' 
#' @param type xxx
#' 
#' @param percent xxx
#' 
#' @param op xxx
#' 
#' @param th xxx
#' 
#' @return An instance of class \code{MSnSet} that have been filtered.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10,]
#' level <- GetTypeDataset(obj)
#' pattern <- 'missing'
#' type <- 'AllCond'
#' percent <- FALSE
#' op <- '=='
#' th <- 2
#' indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, percent, op, th )
#' 
#' @export
#' 
GetIndices_ComplexQueryFiltering <- function(obj, level, pattern, type, percent, op, th ){
  
  if(missing(obj))
    stop ("'obj' is required.")
  if(missing(level))
    stop ("'level' is required.")
  if(missing(pattern))
    stop ("'pattern' is required.")
  if(missing(type))
    tsop ("'type' is required.")
  if(missing(percent))
    stop ("'percent' is required.")
  if(missing(op))
    stop ("'op' is required.")
  if(missing(th))
    stop ("'th' is required.")
  
  
  indices <- NULL
  
  if (!(pattern %in% DAPAR::metacell.def(level)$node && type !='None' && !is.null(type))){
    warning("Either 'pattern' nor 'type' are equal to 'None'")
    return(NULL)
  }
  
  mask <- match.metacell(metadata = Get_qMetadata(obj), 
                         pattern  =pattern, 
                         level = level)
  
  indices <- switch(type,
                    WholeLine = GetIndices_WholeLine(metacell.mask = mask),
                    WholeMatrix = GetIndices_WholeMatrix(metacell.mask = mask,
                                                         op = op, 
                                                         percent = percent, 
                                                         th = th),
                    AllCond = GetIndices_BasedOnConditions(metacell.mask = mask, 
                                                           type = type, 
                                                           conds = Biobase::pData(obj)$Condition, 
                                                           percent = percent, 
                                                           op = op, 
                                                           th = th),
                    AtLeastOneCond = GetIndices_BasedOnConditions(metacell.mask = mask, 
                                                                  type = type,
                                                                  conds =  Biobase::pData(obj)$Condition, 
                                                                  percent = percent,
                                                                  op = op, 
                                                                  th = th)
  )
  
  return(indices)
}




#' @title 
#' Search lines which respects request on one or more conditions.
#' 
#' @description 
#' This function looks for the lines that respect the request in either all conditions
#' or at least one condition.
#' 
#' @param metacell.mask xxx
#' 
#' @param op  String for operator to use. List of operators is available with SymFilteringOperators().
#' 
#' @param percent A boolean to indicate whether the threshold represent an absolute value (percent = FALSE) or
#' a percentage (percent=TRUE).
#' 
#' @param th A floating number which is in the interval [0, 1]
#' 
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing'
#' metacell.mask <- match.metacell(metadata=Get_qMetadata(obj), pattern=pattern, level=level)
#' percent <- FALSE
#' th <- 3
#' op <- '>='
#' ind <- GetIndices_WholeMatrix(metacell.mask, op, percent, th)
#' 
#' @export
#' 
GetIndices_WholeMatrix <- function( metacell.mask, op='==', percent=FALSE, th=0){
  
  #Check parameters
  if (missing(metacell.mask))
    stop("'metacell.mask' is required.")
  if(isTRUE(percent)){
    if (th < 0 || th > 1)
      stop("With percent=TRUE, the threshold 'th' must be in the interval [0, 1].")
  } else {
    th.upbound <- ncol(metacell.mask)
    if (th > th.upbound)
      stop(paste0("Param `th` is not correct. It must be an integer greater than or equal to 0 and less or equal than ",
                  th.upbound))
  }
  
  if (!(op %in% SymFilteringOperators()))
    stop(paste0("'op' must be one of the followinf values: ",
                paste0(SymFilteringOperators(), collapse=' '))
    )
  
  indices <- NULL
  if (isTRUE(percent)) {
    inter <- rowSums(metacell.mask)/ncol(metacell.mask)
    indices <- which(eval(parse(text=paste0("inter", op, th))))
  } else {
    inter <- apply(metacell.mask, 1, sum)
    indices <- which(eval(parse(text=paste0("inter", op, th))))
  }
  
  
  if (length(indices)==0) indices <- NULL
  return(indices)
}


#' @title 
#' Search lines which respects query on all their elements.
#' 
#' @description 
#' This function looks for the lines where each element respect the query.
#' 
#' @param metacell.mask xxx
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[20:30]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing POV'
#' metacell.mask <- match.metacell(metadata=GetMetacell(obj), pattern=pattern, level=level)
#' ind <- GetIndices_WholeLine(metacell.mask)
#'
#'@export
#'
GetIndices_WholeLine <- function(metacell.mask){
  if(missing(metacell.mask))
    stop("'metacell.mask' is missing.")
  
  indices <- which(rowSums(metacell.mask) == ncol(metacell.mask))
  if (length(indices)==0) indices <- NULL
  return (indices)
}


#' @title 
#' Search lines which respects request on one or more conditions.
#' 
#' @description 
#' This function looks for the lines that respect the request in either all conditions
#' or at least one condition.
#' 
#' @param metacell.mask xxx
#' 
#' @param type Available values are:
#' * 'AllCond' (the query is valid in all the conditions),
#' * 'AtLeatOneCond' (the query is valid in at leat one condition.
#' 
#' @param conds xxx
#' 
#' @param percent xxx
#' 
#' @param op  String for operator to use. List of operators is available with SymFilteringOperators().
#' 
#' @param th The theshold to apply
#' 
#' #' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing'
#' metacell.mask <- match.metacell(metadata=Get_qMetadata(obj), pattern=pattern, level=level)
#' type <- 'AllCond'
#' conds <- Biobase::pData(obj)$Condition
#' op <- '>='
#' th <- 2
#' percent <- F
#' ind <- GetIndices_BasedOnConditions(metacell.mask, type, conds, percent, op, th)
#'
#'@export
#'
GetIndices_BasedOnConditions <- function(metacell.mask,
                                         type,
                                         conds, 
                                         percent,
                                         op, 
                                         th){
  
  #Check parameters
  if (missing(metacell.mask))
    stop("'metacell.mask' is missing.")
  if(missing(conds))
    stop("'conds' is missing.")
  if(missing(type))
    stop("'type' is missing.")
  else if (!(type %in% c('AllCond', 'AtLeastOneCond')))
    stop("'type' must be one of the following: AllCond, AtLeastOneCond.")
  if (missing(percent))
    stop("'percent' is missing.")
  if(missing(op))
    stop("'op' is missing.")
  if(missing(th))
    stop("'th' is missing.")
  else if (!(op %in% SymFilteringOperators()))
    stop(paste0("'op' must be one of the followinf values: ",
                paste0(SymFilteringOperators(), collapse=' '))
    )
  
  u_conds <- unique(conds)
  nbCond <- length(u_conds)
  
  
  if(isTRUE(percent)){
    if (th < 0 || th > 1)
      stop("With percent=TRUE, the threshold 'th' must be in the interval [0, 1].")
  } else {
    th.upbound <- min(unlist(lapply(u_conds, function(x) length(which(conds==x)))))
    if (th > th.upbound)
      stop(paste0("Param `th` is not correct. It must be an integer greater than or equal to 0 and less or equal than ",
                  th.upbound))
  }
  
  
  
  
  indices <- NULL
  s <- matrix(rep(0, nrow(metacell.mask)*nbCond),
              nrow=nrow(metacell.mask),
              ncol=nbCond)
  
  for (c in 1:nbCond) {
    ind.cond <- which(conds == u_conds[c])
    inter <- rowSums(metacell.mask[, ind.cond])
    if (isTRUE(percent))
      inter <- inter/length(ind.cond)
    s[,c] <- eval(parse(text=paste0("inter", op, th)))
  }
  
  indices <- switch(type,
                    AllCond = which(rowSums(s) == nbCond),
                    AtLeastOneCond = which(rowSums(s) >= 1)
  )
  
  return(indices)
}




#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @description Returns the \code{SummarizedExperiment} object with a extra column in \code{rowData()}.
#' The extra column indicates feature(s) to delete in 1 and 0 otherwise.
#' The user chooses the threshold, either the percentage of NAS per row (default)
#' or the number of samples containing NAs per row allowed. Then, 
#' the filter tags the lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#'
#' The different methods are :
#' - "WholeMatrix": given a threshold \code{th}, only the lines that contain
#'   at least \code{th} values are kept.
#' - "AllCond": given a threshold \code{th}, only the lines which contain
#'   at least \code{th} values for each of the conditions are kept.
#' - "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#'   at least \code{th} values, and for at least one condition, are kept.
#'
#' @param object An object of class \code{QFeatures}
#' 
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' 
#' @param th Either a numeric between 0 and 1 where only the lines which contain
#' at least \code{th}% of non NA values are kept.
#' Or a integer between 0 and maximum number of samples for 'WholeMatrix' and 
#' between 0 and maximum number of replicate for "AllCond" and "AtLeastOneCond", where
#' only the lines which contain at least \code{th} values are kept.
#' 
#' @param percent TRUE by default. When FALSE, use the number of samples
#' 
#' @return The object of class \code{SummarizedExperiment} with extra column in rowData
#' indicating 1 for the lines to remove, else 0.
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept
#' res <- MVrowsTagToOne(object, type = "AtLeastOneCond", th=0.7, percent=TRUE)
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' 
MVrowsTagToOne <- function(object, type, th=0, percent=TRUE) {
  
  if (is.null(object)) { return(NULL) }
  
  if (missing(type))
    stop("'type' is required")
  
  if (!(type %in% c('None', 'EmptyLines', 'WholeMatrix', 'AtLeastOneCond', 'AllCond')))
    stop("'type' is not one of: 'None', 'EmptyLines', 'WholeMatrix', 'AtLeastOneCond', 'AllCond'")
  
  newColName <- "tagNA"
  i <- length(experiments(object))
  
  ## Create a fictive column for each assays otherwise the filter method
  ## of QFeatures will truncate the object
  for (k in 1:length(experiments(object)))
    rowData(object[[k]]) <- setNames(cbind(rowData(object[[k]]),tmp=rep(0, nrow(object[[k]]))),
                                     c(names(rowData(object[[k]])),newColName)
    )
  
  
  sampleTab <- colData(object)
  
  
  # Filtration
  keepThat <- df <- NULL
  df <- as.data.frame(assay(object[[i]]))
  
  
  if (type == "None") {
    keepThat <- seq(1:nrow(df))
  } else if (type == "EmptyLines") {
    keepThat <- which(apply(!is.na(df), 1, sum) >= 1)
  } else if (type == "WholeMatrix") {
    if (isTRUE(percent)) {
      keepThat <- which(rowSums(!is.na(df))/ncol(df) >= th) 
    } else {
      keepThat <- which(apply(!is.na(df), 1, sum) >= th)
    }
  } else if (type == "AtLeastOneCond" || type == "AllCond") {
    
    if (is.null(sampleTab)) { return(NULL) }
    
    conditions <- unique(sampleTab[['Condition']])
    nbCond <- length(conditions)
    s <- matrix(rep(0, nrow(df)*nbCond),nrow=nrow(df), ncol=nbCond)
    
    if (isTRUE(percent)) {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        s[,c] <- (rowSums(!is.na(df[,ind]))/length(ind)) >= th
      }
    } else {
      for (c in 1:nbCond) {
        ind <- which(sampleTab[['Condition']] == conditions[c])
        if (length(ind) == 1){
          s[,c] <- !is.na(df[,ind]) >= th 
        }
        else {
          s[,c] <- (apply(!is.na(df[,ind]), 1, sum)) >= th
        }
      }
    }
    
    if (type == "AllCond")
      keepThat <- which(rowSums(s) == nbCond)
    else if (type == "AtLeastOneCond")
      keepThat <- which(rowSums(s) >= 1)
  }
  
  
  rowData(object[[i]])[-keepThat, newColName] <- 1 
  
  return(object)
}



# TODO get ideas in this function to update the previous one. This function is in DAPAR
#' 
#' #' Returns the indices of the lines of \code{exprs()} table to delete w.r.t. 
#' #' the conditions on the number of missing values.
#' #' The user chooses the minimum amount of intensities that is acceptable and
#' #' the filter delete lines that do not respect this condition.
#' #' The condition may be on the whole line or condition by condition.
#' #' 
#' #' The different methods are :
#' #' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' #' at least \code{th} values are kept.
#' #' "AllCond": given a threshold \code{th}, only the lines which contain
#' #' at least \code{th} values for each of the conditions are kept.
#' #' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' #' at least \code{th} values, and for at least one condition, are kept.
#' #' 
#' #' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' #' 
#' #' @param obj An object of class \code{MSnSet} containing
#' #' quantitative data.
#' #' 
#' #' @param percent TRUE or FALSE. Default is FALSE..
#' #' 
#' #' @param condition Method used to choose the lines to delete.
#' #' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' #' 
#' #' @param threshold An integer value of the threshold if percent is FALSE. Otherwise, a floating
#' #' number between 0 and 1.
#' #' 
#' #' @return An vector of indices that correspond to the lines to keep.
#' #' 
#' #' @author Enora Fremy, Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' mvFilterGetIndices(Exp1_R25_pept, condition = "WholeMatrix", threshold=2)
#' #' mvFilterGetIndices(Exp1_R25_pept, condition = "EmptyLines")
#' #' mvFilterGetIndices(Exp1_R25_pept, condition = "WholeMatrix", percent=TRUE, threshold=0.5)
#' #' 
#' #' @export
#' #' 
#' mvFilterGetIndices <- function(obj,
#'                                percent = FALSE,
#'                                condition = 'WholeMatrix', 
#'                                threshold = NULL){
#'   #Check parameters
#'   paramtype<-c("None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond")
#'   if (!(condition %in% paramtype)){
#'     warning("Param `type` is not correct.")
#'     return (NULL)
#'   }
#'   
#'   if (condition != 'EmptyLines')
#'     if (!(percent %in% c(T, F))){
#'       warning("Param `type` is not correct.")
#'       return (NULL)
#'     } else {
#'       if (!isTRUE(percent)){
#'         paramth <- c(seq(0, nrow(Biobase::pData(obj)), 1))
#'         if (!(threshold %in% paramth)){
#'           warning(paste0("Param `threshold` is not correct. It must an integer greater than or equal to 0 and less or equal than ",
#'                          nrow(Biobase::pData(obj))))
#'           return (NULL)
#'         }
#'       } else {
#'         if (threshold < 0 || threshold > 1){
#'           warning("Param `threshold` is not correct. It must be greater than 0 and less than 1.")
#'           return (NULL)
#'         }
#'       }
#'     }
#'   
#'   keepThat <- NULL
#'   if (is.null(obj@experimentData@other$OriginOfValues)){
#'     data <- Biobase::exprs(obj)
#'     warning('The dataset contains no slot OriginOfValues in which to search for indices. The search will
#'             be proceeded in the intensities tab based on NA values')
#'   } else {
#'     data <- dplyr::select(Biobase::fData(obj),
#'                           obj@experimentData@other$OriginOfValues)
#'   }
#'   
#'   if (condition == "None") {
#'     keepThat <- seq(1:nrow(data))
#'   } else if (condition == "EmptyLines") {
#'     keepThat <- which(apply(!DAPAR::is.MV(data), 1, sum) >= 1)
#'   } else if (condition == "WholeMatrix") {
#'     if (isTRUE(percent)) {
#'       keepThat <- which(rowSums(!DAPAR::is.MV(data))/ncol(data) >= threshold) 
#'     } else {
#'       keepThat <- which(apply(!DAPAR::is.MV(data), 1, sum) >= threshold)
#'     }
#'   } else if (condition == "AtLeastOneCond" || condition == "AllCond") {
#'     
#'     conditions <- unique(Biobase::pData(obj)$Condition)
#'     nbCond <- length(conditions)
#'     keepThat <- NULL
#'     s <- matrix(rep(0, nrow(data)*nbCond),
#'                 nrow=nrow(data),
#'                 ncol=nbCond)
#'     
#'     if (isTRUE(percent)) {
#'       for (c in 1:nbCond) {
#'         ind <- which(Biobase::pData(obj)$Condition == conditions[c])
#'         s[,c] <- (rowSums(!DAPAR::is.MV(data[,ind]))/length(ind)) >= threshold
#'       }
#'     } else {
#'       for (c in 1:nbCond) {
#'         ind <- which(Biobase::pData(obj)$Condition == conditions[c])
#'         if (length(ind) == 1){
#'           s[,c] <- (!DAPAR::is.MV(data[,ind]) >= threshold) 
#'         }
#'         else {
#'           s[,c] <- (apply(!DAPAR::is.MV(data[,ind]), 1, sum)) >= threshold
#'         }
#'       }
#'     }
#'     
#'     switch(condition,
#'            AllCond = keepThat <- which(rowSums(s) == nbCond),
#'            AtLeastOneCond = keepThat <- which(rowSums(s) >= 1)
#'     )
#'   }
#'   
#'   return(keepThat)
#' }
#' 



#' @title
#' Filter missing values by proportion
#'
#' @description Remove lines in the data according to the proportion of missing
#' values. This proportion is calculated differently depending on whether we
#' want a certain proportion of missing values (NA) to remain on:
#' * the entire matrix, regardless of the conditions: the rows containing a
#' proportion of NA equal or below the threshold will be kept.
#' * all the conditions: the lines for which all the conditions have a NA
#' proportion equal to or less than the fixed proportion will be kept.
#' * at least one condition: the lines for which at least one condition is
#' equal to or less than the fixed proportion of NA will be kept.
#'
#' @param obj  An object of class \code{SummarizedExperiment} containing quantitative data
#' and phenotype data.
#' 
#' @param sTab xxxx
#' 
#' @param int.prop float between 0 and 1 corresponding to the proportion
#' of intensities to keep in the lines.
#' 
#' @param mode character string. Four possibilities corresponding to the
#' description above: "None", WholeMatrix", "AllCond" and "AtLeastOneCond".
#' 
#' @return the object given as input but with the lines not respecting the
#' proportion of NA requested in less.
#' 
#' @author H?l?ne Borges, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_prot, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' sTab <- colData(obj)
#' 
#' obj <- MVrowsTagToOne_HB(obj, sTab, mode = 'AtLeastOneCond', int.prop = 0.7)
#' na_filter_percent <- VariableFilter(field = "tagNA", value = "0", condition = "==")
#' obj <- filterFeaturesSam(object = obj, i = 2, name = 'na_filter_percent', filter=na_filter_percent)
#' obj <- removeAdditionalCol(obj, "tagNA")
#' 
#' @export
#' 
#' @importFrom stringr str_glue
#' @import dplyr
#' @importFrom tidyr pivot_longer %>%
#' @importFrom methods is
#' 
MVrowsTagToOne_HB <- function(obj, sTab, int.prop, mode = "None"){
  
  if(class(obj) != 'QFeatures')
    stop("'obj' must be of class 'SummarizedExperiment'")
  
  if(missing(sTab))
    stop("'sTab' is required")
  
  # check if mode is valid
  if(!(mode %in% c("None","WholeMatrix", "AllCond", "AtLeastOneCond"))){
    stop(stringr::str_glue("Wrong mode: {mode} is not a valid string.
                     Please choose between 'None', 'WholeMatrix', 'AllCond' or 'AtLeastOneCond'.",
                           call. =FALSE))
  }
  # check if int.prop is valid
  if(missing(int.prop))
    stop("'int.prop' is needed")
  else 
  {
    if(!methods::is(int.prop, "numeric" )){
      stop(stringr::str_glue("Wrong parameter: int.prop needs to be numeric"))
    } else if(!dplyr::between(int.prop, 0, 1)){
      stop(stringr::str_glue("Wrong parameter: int.prop must be between 0 and 1"))
    }
  }
  
  
  newColName <- "tagNA"
  i <- length(experiments(obj))
  
  ## Create a fictive column for each assays otherwise the filter method
  ## of QFeatures will truncate the object
  for (k in 1:length(experiments(obj)))
    rowData(obj[[k]]) <- setNames(cbind(rowData(obj[[k]]),tmp=rep(0, nrow(obj[[k]]))),
                                  c(names(rowData(obj[[k]])),newColName)
    )
  
  
  
  
  
  
  
  print(stringr::str_glue("Choosen proportion of intensities to be present: {int.prop}"))
  print(stringr::str_glue("Choosen mode: {mode}"))
  intensities <- assay(obj)
  sTab$Condition <- as.factor(sTab$Condition)
  
  intensities_t <- as.data.frame(t(intensities))
  intensities_t <- dplyr::bind_cols(intensities_t,
                                    condition = sTab$Condition,
                                    sample = rownames(intensities_t))
  tbl_intensities <- dplyr::as_tibble(intensities_t, rownames = NA)
  longer_intensities <- tbl_intensities %>%
    tidyr::pivot_longer(-c(condition, sample), names_to = "feature", values_to = "intensity")
  # group_by does not keep the initial order when it is not a factor so to keep
  # the protein order, we cheat by transforming feature into a factor.
  longer_intensities$feature <- factor(longer_intensities$feature,
                                       levels = unique(longer_intensities$feature))
  if(mode == "None"){
    to_keep <- 1:nrow(obj)
  }else if(mode == "WholeMatrix"){
    nb_samples <- ncol(intensities)
    threshold <- ceiling(nb_samples*int.prop)
    print(stringr::str_glue("missing value threshold {threshold}"))
    # for each feature (protein / peptide) we count the number of intensities present
    feat_grp <- longer_intensities %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(non_na = sum(!is.na(intensity)))
    to_keep <- which(feat_grp$non_na >= threshold)
    
  }else if(mode == "AllCond" || mode == "AtLeastOneCond"){
    workforces <- longer_intensities %>%
      dplyr::group_by(feature, condition) %>%
      dplyr::count(condition)
    # the number of samples per condition
    workforces <- workforces$n[seq_len(length(levels(sTab$Condition)))]
    
    # for each condition of each feature, we count the number of intensities present
    feat_grp <- longer_intensities %>%
      dplyr::group_by(feature, condition) %>%
      dplyr::summarise(non_na = sum(!is.na(intensity)))
    # the threshold for each condition
    thresholds <- ceiling(workforces*int.prop)
    print(stringr::str_glue("for condition {unique(levels(longer_intensities$condition))} number of samples is {workforces}, so missing value threshold is {thresholds} "))
    # For each feature, each condition is compared with its respective
    # threshold, we put 0 if the protein has a number of intensities lower than
    # the threshold of the corresponding condition, and 1 otherwise
    check_th <- feat_grp %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(non_na = dplyr::case_when(
        non_na < thresholds ~ 0,
        TRUE ~ 1
      )) %>%
      dplyr::ungroup()
    # if it is AllCond then we must find the features for which all the conditions
    # respect the threshold
    if(mode == "AllCond"){
      all_cond_ok <- check_th %>%
        dplyr::group_by(feature) %>%
        dplyr::filter(all(non_na ==1)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
      all_cond_ok$feature <- as.character(all_cond_ok$feature)
      to_keep <- which(rownames(obj[[length(experiments(obj))]]) %in% all_cond_ok$feature)
    }else if(mode == "AtLeastOneCond"){
      # if it is AtLeastOneCond then we must find the features for which at
      # least one condition that respects the threshold
      any_cond_ok <- check_th %>%
        dplyr::group_by(feature) %>%
        dplyr::filter(any(non_na ==1)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
      any_cond_ok$feature <- as.character(any_cond_ok$feature)
      to_keep <- which(rownames(obj[[length(experiments(obj))]]) %in% any_cond_ok$feature)
    }
  }
  print(stringr::str_glue("There were initially {nrow(intensities)} features.
                 After filtering out the missing values, {nrow(to_keep)} remain."))
  
  
  rowData(obj[[length(experiments(obj))]])[-to_keep, newColName] <- 1 
  
  return(obj)
}











#' @title Restore the rowData() before MVrowsTagToOne
#' 
#' @param object An object of class \code{QFeatures}
#' 
#' @param colToRemove Column to remove
#' 
#' @return The \code{SummarizedExperiment} object without the rowData colToRemove column
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[[2]]
#' sampleTab <- colData(Exp1_R25_pept)
#' obj <- MVrowsTagToOne(obj=object, i=2, sampleTab, newColName="LinesKept", type="AllCond", th=2)
#' obj <- removeAdditionalCol(obj, i=2, colToRemove="LinesKept" )
#' 
#' @export
#' 
#' @import SummarizedExperiment
#' 
removeAdditionalCol <- function(object, colToRemove=NULL) {
  
  if(is.null(object)) { return(NULL) }
  if(is.null(colToRemove)){return(NULL)} 
  colToRemove <- "tagNA"
  
  for (k in 1:length(experiments(object)))
  {
    if( is.na(match(colToRemove,names(rowData(object[[k]])))) ) {
      print(paste0("Warning: ",colToRemove," isn't a column name of rowData(obj)"))
    }
    else 
      rowData(object[[k]]) <- rowData(object[[k]])[,-which(colnames(rowData(object[[k]]))==colToRemove)]
  }
  
  return(object)
  
}








#' getPourcentageOfMV <- function(obj)
#' Voir nNA(object, i) de la classe QFeatures






#' @title Barplot of proportion of contaminants and reverse
#' 
#' @param lDataset The total length (number of rows) of the dataset
#' 
#' @param nBoth The number of both contaminants and reverse identified in the dataset.
#' 
#' @param nCont The number of contaminants identified in the dataset.
#' 
#' @param nRev The number of reverse entities identified in the dataset.
#' 
#' @return A barplot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' proportionConRev_HC(100, 10, 20)
#' 
#' @export
#' 
#' @import highcharter
#' 
proportionConRev_HC <- function(lDataset, nBoth = 0, nCont=0, nRev=0){
  
  if(missing(lDataset))
    stop("'lDataset' is missing.")
  
  
  total <- nBoth + nCont + nRev + lDataset
  pctGood <- 100 * round(lDataset/total,  digits=4)
  pctBoth <- 100 * round(nBoth/total,  digits=4)
  pctContaminants <- 100 * round(nCont/total,  digits=4)
  pctReverse <- 100 * round(nRev/total,  digits=4)
  
  counts <- c(lDataset, nCont, nRev, nBoth)
  slices <- c(pctGood, pctContaminants, pctReverse ,pctBoth) 
  lbls <- c("Quantitative data", "Contaminants", "Reverse", "Both contaminants & Reverse")
  lbls <- paste(lbls, " (", counts, " lines)", sep="") 
  
  mydata <- data.frame(test=c(pctGood, pctContaminants, pctReverse ,pctBoth))
  
  highchart() %>% 
    dapar_hc_chart(chartType = "bar") %>% 
    hc_yAxis(title = list(text = "Pourcentage")) %>% 
    hc_xAxis(categories=lbls) %>% 
    hc_legend(enabled = FALSE) %>%
    hc_plotOptions(column = list(
      dataLabels = list(enabled = TRUE),
      stacking = "normal",
      enableMouseTracking = FALSE)
    ) %>% 
    hc_add_series(data  = mydata$test,
                  dataLabels = list(enabled = TRUE, format='{point.y}%'),
                  colorByPoint = TRUE) %>%
    dapar_hc_ExportMenu(filename = "contaminants")
  
  
}