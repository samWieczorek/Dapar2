

#' @title Quantitative cell metadata scopes for filtering
#'
#' @export
#' 
#' @rdname qMetadata-filter
#'
#' @return NA
#'
qMetadataFilteringScope <- function()
  c("None", "WholeLine", "WholeMatrix", "AllCond", "AtLeastOneCond")



#' @title Operators for complex queries
#'
#' @export
#' 
#' @rdname qMetadata-filter
#'
#' @return NA
#'
SymFilteringOperators <- function()
  c('<=','<', '>=', '>', '==', '!=')



#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either all conditions
#' or at least one condition.
#'
#' @param mask xxx
#'
#' @param op  String for operator to use. List of operators is available with SymFilteringOperators().
#'
#' @param percent A boolean to indicate whether the threshold represent an absolute value (percent = FALSE) or
#' a percentage (percent=TRUE).
#'
#' @param th A floating number which is in the interval [0, 1]
#'
#' @return NA
#'
#' @examples
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing'
#' mask <- match.qMetadata(metadata=Get_qMetadata(obj), pattern=pattern, level=level)
#' percent <- FALSE
#' th <- 3
#' op <- '>='
#' ind <- GetIndices_WholeMatrix(mask, op, percent, th)
#'
#' @export
#' 
#' @rdname qMetadata-filter
#'
qMetadatWholeMatrix <- function(object, cmd, pattern, percent, th, operator){
  if(missing(object))
    return(NULL)
  if(missing(cmd))
    return(object)
  if(missing(pattern))
    return(object)
  if(missing(percent))
    return(object)
  if(missing(th))
    return(object)
  if(missing(operator))
    return(object)
  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('qMetadata' %in% names(rowData(object)))
  
  stopifnot(!is.null(cmd) || (cmd %in% c('keep', 'delete')))
  
  indices <- NULL
  level <- typeDataset(object)
  
  
  #Check parameters
  mask <- match.qMetadata(df = qMetadata(object),
                          pattern = pattern,
                          level = level)
  if(isTRUE(percent)){
    if (th < 0 || th > 1)
      stop("With percent=TRUE, the threshold 'th' must be in the interval [0, 1].")
  } else {
    if (th > ncol(mask))
      stop(paste0("Param `th` is not correct. It must be an integer greater than or equal to 0 and less or equal than ",
                  ncol(mask)))
  }
  
  if (!(operator %in% SymFilteringOperators()))
    stop(paste0("'op' must be one of the followinf values: ",
                paste0(SymFilteringOperators(), collapse=' '))
    )
  
  
  
  indices <- NULL
  if (isTRUE(percent)) {
    inter <- rowSums(mask)/ncol(mask)
    indices <- which(eval(parse(text = paste0("inter", operator, th))))
  } else {
    inter <- apply(mask, 1, sum)
    indices <- which(eval(parse(text = paste0("inter", operator, th))))
  }
  
  if (length(indices) >0)
    object <- switch(cmd, 
                     keep = object[indices,],
                     delete = object[-indices,])
  
  return (object)
}


#' @title
#' Search lines which respects query on all their elements.
#'
#' @description
#' This function looks for the lines where each element respect the query.
#'
#' @param mask xxx
#'
#' @return NA
#'
#' @examples
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[20:30]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing POV'
#' mask <- match.qMetadata(metadata=GetqMetadata(obj), pattern=pattern, level=level)
#' ind <- GetIndices_WholeLine(mask)
#'
#'@export
#'
#' @rdname qMetadata-filter
#'
qMetadataWholeLine <- function(object, cmd, pattern){
  if(missing(object))
    return(NULL)
  if(missing(cmd))
    return(object)
  if(missing(pattern))
    return(object)
  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  stopifnot('qMetadata' %in% names(rowData(object)))
  stopifnot(!is.null(cmd) || (cmd %in% c('keep', 'delete')))
  
  indices <- NULL
  level <- typeDataset(object)
  
  if (!(pattern %in% qMetadata.def(level)$node)){
    warning("Either 'pattern' nor 'type' are equal to 'None'")
    return(NULL)
  }
  
  mask <- match.qMetadata(df = qMetadata(object),
                          pattern = pattern,
                          level = level)
  
  
  indices <- which(rowSums(mask) == ncol(mask))
  if (length(indices) >0)
    object <- switch(cmd, 
                     keep = object[indices,],
                     delete = object[-indices,])
  
  return (object)
}


#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either all conditions
#' or at least one condition.
#'
#' @param mask xxx
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
#' @return NA
#'
#' #' @examples
#' Exp1_R25_pept <- readRDS(system.file("data", 'Exp1_R25_pept.rda', package="DaparToolshedData"))
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- obj@experimentData@other$typeOfData
#' pattern <- 'missing'
#' mask <- match.qMetadata(metadata=qMetadata(obj), pattern=pattern, level=level)
#' type <- 'AllCond'
#' conds <- Biobase::pData(obj)$Condition
#' op <- '>='
#' th <- 2
#' percent <- FALSE
#' ind <- GetIndices_BasedOnConditions(mask, type, conds, percent, op, th)
#'
#' @export
#'
#' @rdname qMetadata-filter
#'
qMetadataOnConditions <- function(object,
                                  cmd,
                                  mode,
                                  pattern,
                                  conds,
                                  percent,
                                  operator,
                                  th){
  
  #Check parameters
  if (missing(pattern))
    stop("'pattern' is missing.")
  if(missing(conds))
    stop("'conds' is missing.")
  if(missing(mode))
    stop("'type' is missing.")
  else if (!(mode %in% c('AllCond', 'AtLeastOneCond')))
    stop("'mode' must be one of the following: AllCond, AtLeastOneCond.")
  if (missing(percent))
    stop("'percent' is missing.")
  
  if(missing(th))
    stop("'th' is missing.")
  else
    stopifnot(th >= 0 && th <= 1)
  
  if(missing(operator))
    stop("'operator' is missing.")
  else if (!(operator %in% SymFilteringOperators()))
    stop(paste0("'operator' must be one of the followinf values: ",
                paste0(SymFilteringOperators(), collapse=' '))
    )
  
  u_conds <- unique(conds)
  nbCond <- length(u_conds)
  
  if(!isTRUE(percent)){
    th.upbound <- min(unlist(lapply(u_conds, 
                                    function(x) length(which(conds==x)))))
    if (th > th.upbound)
      stop(paste0("Param `th` is not correct. It must be an integer greater than or equal to 0 and less or equal than ",
                  th.upbound))
  }
  
  
  mask <- match.qMetadata(df = qMetadata(object),
                          pattern = pattern,
                          level = typeDataset(object)
  )
  
  indices <- NULL
  s <- matrix(rep(0, nrow(mask)*nbCond),
              nrow = nrow(mask),
              ncol = nbCond)
  
  for (c in seq_len(nbCond)) {
    ind.cond <- which(conds == u_conds[c])
    inter <- rowSums(mask[, ind.cond])
    if (isTRUE(percent))
      inter <- inter/length(ind.cond)
    s[,c] <- eval(parse(text=paste0("inter", operator, th)))
  }
  
  indices <- switch(mode,
                    AllCond = which(rowSums(s) == nbCond),
                    AtLeastOneCond = which(rowSums(s) >= 1)
  )
  
  if (length(indices) >0)
    object <- switch(cmd, 
                     keep = object[indices,],
                     delete = object[-indices,])
  
  return (object)
}
