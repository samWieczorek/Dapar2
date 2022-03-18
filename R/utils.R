
#' @title Standardize names
#'
#' @description 
#' 
#' Replace ".", ' ', '-' in `character()` by '_' to be compliant
#' with functions of [Shinyjs], [Shiny]
#'
#' @param x A `character()` to be processed
#'
#' @return A `character()` of the same length as 'x' with modified
#' names.
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples
#' 
#' ReplaceSpecialChars(c('foo.1', 'foo-2', 'foo 3'))
#'
ReplaceSpecialChars <- function(x){
  if (is.null(x))
    return(x)
  
  for (char in c('.', ' ', '-'))
    x <- gsub(char, '_', x, fixed=TRUE)
  
  x
}


#' @title Version number of Prostar suite
#' 
#' @description 
#' 
#' This function gives the version number of the packages of the 
#' Prostar suite and which propose data processing. This information
#' can be useful if the user wants to publish its works or to
#' rerun a data processing pipeline tin a given set of conditions. 
#' The packages which are concerned are [Prostar], [DaparToolshed] 
#' and [DaparToolshedData]
#' 
#' @return A `list(3)`
#' 
#' @rdname ProstarVersions
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' ProstarVersions()
#' 
#' @export
#' 
ProstarVersions <- function(){
  
  Prostar <- DaparToolshed <- DaparToolshedData <- Magellan <- NA
  
  tryCatch({
    find.package("Prostar")
    Prostar <- Biobase::package.version('Prostar')
  },
  error = function(e) Prostar <- NA
  )
  
  tryCatch({
    find.package("DaparToolshed")
    DaparToolshed <- Biobase::package.version('DaparToolshed')
  },
  error = function(e) DaparToolshed <- NA
  )
  
  tryCatch({
    find.package("DaparToolshedData")
    DaparToolshed <- Biobase::package.version('DaparToolshedData')
  },
  error = function(e) DaparToolshedData <- NA
  )

  
  list(Prostar = Prostar,
       DaparToolshed = DaparToolshed,
       DaparToolshedData = DaparToolshedData)
}










#' @title Number of empty lines in the data
#' 
#' @description 
#' 
#' This function counts the number of empty lines (all elements
#' are equal to NA).
#' 
#' @param df A `data.frame`.
#' 
#' @return A `integer(1)`
#' 
#' @export
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' data(ft)
#' nEmptyLines(assay(ft, 1))
#' 
nEmptyLines <- function(df)
  sum(apply(is.na(as.matrix(df)), 1, all))






#' @title Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#' 
#' @param data A data.frame
#' 
#' @param type The value to search in the dataframe
#' 
#' @return A boolean dataframe 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' data(ft)
#' obj <- ft[[1]]
#' data <- qMetadata(obj)
#' is.OfType(as.data.frame(data), "MEC")
#' 
#' @export
#' 
is.OfType <- function(data, type){
  stopifnot(inherits(data, 'data.frame'))
  return (type == data)
}







#' @title Returns the possible number of values in lines in the data
#' 
#' @param obj An object of class \code{QFeatures}
#' 
#' @param i The indice of the dataset (SummarizedExperiment) in the object
#' 
#' @param type WholeMatrix, AllCond or AtLeastOneCond
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek, Enora Fremy
#' 
#' @examples
#' data(ft)
#' res <- getListNbValuesInLines(ft, 1)
#' 
#' @export
#' 
#' @importFrom S4Vectors sort
#' @importFrom SummarizedExperiment colData
#' 
getListNbValuesInLines <- function(obj, i, type="WholeMatrix"){
  
  if (is.null(obj)){return(NULL)}
  if(missing(i))
    stop("'i' is required")
  if (!(type %in% c('None', 'WholeMatrix', 'AtLeastOneCond', 'AllCond')))
    stop("'type' is not one of: 'None', 'WholeMatrix', 'AtLeastOneCond', 'AllCond'")
  
  data <- as.data.frame(SummarizedExperiment::assay(obj[[i]]))
  conds <- colData(obj)[['Condition']]
  switch(type,
         WholeMatrix= {
           ll <- unique(ncol(data) - apply(is.na(data), 1, sum))
         },
         AllCond = {
           tmp <- NULL
           for (cond in unique(conds)){
             tmp <- c(tmp, length(which(conds == cond)))
           }
           ll <- seq(0,min(tmp))
         },
         AtLeastOneCond = {
           tmp <- NULL
           for (cond in unique()){
             tmp <- c(tmp, length(which(conds == cond)))
           }
           ll <- seq(0,max(tmp))
         }
  )
  
  return (sort(ll))
}







#' @title Retrieve the indices of non-zero elements in sparse matrices
#' 
#' @description 
#' 
#' This function retrieves the indices of non-zero elements in sparse matrices
#' of class dgCMatrix from package Matrix. This function is largely inspired from 
#' the package \code{RINGO}
#' 
#' @param x A sparse matrix of class dgCMatrix
#' 
#' @return A two-column matrix
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(Matrix)
#' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow = TRUE, sparse=TRUE)
#' res <- nonzero(mat)
#' 
#' @export
#' 
nonzero <- function(x){
  ## function to get a two-column matrix containing the indices of the
  ### non-zero elements in a "dgCMatrix" class matrix
  
  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0))
    return(matrix(0, nrow=0, ncol=2,
                  dimnames=list(character(0), c("row","col"))))
  res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  colnames(res) <- c("row", "col")
  res <- res[x@x != 0, , drop = FALSE]
  return(res)
}



#' @param fname xxx
#' @export
GetExtension <- function(fname)
  strsplit(fname, '.', TRUE)[[1]][2]