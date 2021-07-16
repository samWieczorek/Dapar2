#' @title List normalization methods
#' 
#' @param withTracking A boolean to indicate methods with tracking feature
#' 
#' @value A list of methods
#' 
#' @name normalizeMethods.dapar
#' 
#' @export
#'
normalizeMethods.dapar <- function(withTracking=FALSE){
  if (isTRUE(withTracking))
    c("SumByColumns", 
      "QuantileCentering", 
      "MeanCentering")
else
  c("GlobalQuantileAlignment", 
    "SumByColumns",
    "QuantileCentering", 
    "MeanCentering",
    "LOESS", 
    "vsn")
}



#' @title Check the validity of the experimental design.
#'
#' @description This manual page describes the computation of statistical test using [QFeatures] objects. In the following
#' functions, if `object` is of class `QFeatures`, and optional assay
#' index or name `i` can be specified to define the assay (by name of
#' index) on which to operate.
#'
#' The following functions are currently available:
#'
#' - `GlobalQuantileAlignment` xxx
#' - `SumByColumns` xxx
#' - `QuantileCentering` xxx
#' - `MeanCentering` xxxx
#' - `LOESS` xxx
#' - `vsn` xxx
#' 
#'
#' @details
#' #' Provides several methods to normalize quantitative data from a numeric matrix.
#' They are organized in six main families : GlobalQuantileAlignement, sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
#' For the first family, there is no type.
#' For the five other families, two type categories are available : 
#' "Overall" which means that the value for each protein (ie row) is computed over all the samples ;
#' "within conditions" which means that the value for each protein (ie row) is computed condition by condition.
#' Normalization can be based only on a subset of protein.
#' Proteins with NA in this subset are ignored.
#'
#' 
#' @return An instance of class 'SummerizedExperiment' where the quantitative
#' data in the \code{array()} tab has been normalized.
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept

#' conds <- colData(object)$Condition
#' object <- normalizeD(object, i=2, name = "norm", method='SumByColumns', conds=conds, type='overall')
#' 
"normalizeD"


#' @param  object An object of class `SummarizedExperiment`.
#' 
#' @param method xxxxx. See xxx for methods available.
#' 
#' @param withTracking xxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname normalizeD
#' 
setMethod("normalizeD", "SummarizedExperiment",
          function(object,
                   method,
                   withTracking=FALSE,
                   ...) {
            
            if(!(method %in% normalizeMethods.dapar(withTracking)))
              stop("'method' is not correct. See 'normalizeMethods.dapar()' for details.")
            
            e <- do.call(method, list(assay(object), ...))
            
            rownames(e) <- rownames(assay(object))
            colnames(e) <- colnames(assay(object))
            assay(object) <- e
            object
          })



#' @param  object An object of class `QFeatures`.
#' 
#' @param method xxxxx. See xxx for methods available.
#' 
#' @param i xxxx
#' 
#' @param name xxxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @importFrom QFeatures addAssay addAssayLinkOneToOne
#' 
#' @rdname normalizeD
#' 
setMethod("normalizeD", "QFeatures",
          function(object, i, name = "normalizedAssay", ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            
            argg <- c(as.list(environment()), list(...))
            
            tmp <-  normalizeD(object[[i]], ...)
            metadata(tmp)$Params <- argg[-match(c('object', 'i', 'name'), names(argg))]
            
            object <- QFeatures::addAssay(object,
                                         tmp,
                                         name)
            QFeatures::addAssayLinkOneToOne(object, from = i, to = name)
          })



#' @title Normalisation GlobalQuantileAlignement
#' 
#' @param qData xxxx
#' 
#' @return A normalized numeric matrix
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' normalized <- GlobalQuantileAlignment(qData)
#' 
#' @export
#' 
#' @importFrom preprocessCore normalize.quantiles
#' 
GlobalQuantileAlignment <- function(qData) {
  e <- preprocessCore::normalize.quantiles(as.matrix(qData))
  return(e)
}


#' @title Normalisation SumByColumns
#' 
#' @param qData xxxx
#' 
#' @param conds xxx
#' 
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' 
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' 
#' @return A normalized numeric matrix
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' library(QFeatures)
#' qData <- assay(Exp1_R25_pept[['original_log']])
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' normalized <- SumByColumns(qData, conds, type="within conditions", subset.norm=1:10)
#' 
#' @export
#' 
#' @importFrom stats median
#' 
SumByColumns <- function(qData, 
                         conds = NULL, 
                         type = NULL, 
                         subset.norm = NULL) {

  if (missing(qData))
    stop("'qData' is missing.")

if( missing(conds) || is.null(conds))
    stop("'conds' is required")
  
  if( missing(type) || is.null(type)){
    stop("'type' is required")
    return()
  }
  
  
  if (!(type %in% c('overall', 'within conditions')))
    stop("'type' must have one of the following values: 'overall', 'within conditions'")
  
  
  qData <- as.matrix(qData)
  
  e <- 2^qData
  
  if(is.null(subset.norm) || length(subset.norm) < 1){
    subset.norm = 1:nrow(qData)
  }
  
  if (type == "overall"){
    
    
      if(length(subset.norm) == 1)
        sum_cols = e[subset.norm,]
      else
        sum_cols <- colSums(e[subset.norm,], na.rm=TRUE)
      
    
    for ( i in 1:nrow(e))
      e[i, ] <- (e[i, ] / sum_cols)*(stats::median(sum_cols))
  } else if (type == "within conditions"){
    
    for (l in unique(conds)) {
      indices <- which(conds== l)
      
          if(length(subset.norm)==1)
            sum_cols=e[subset.norm,indices]
          else
            sum_cols <- colSums(e[subset.norm,indices], na.rm=TRUE)
      
      for (i in 1:nrow(e))
        e[i,indices] <- (e[i,indices]/sum_cols) * stats::median(sum_cols)
    }
  }
  e <- log2(e)
  return(e)
}


#' @title Normalisation QuantileCentering
#' 
#' @param qData xxx
#' 
#' @param conds xxx
#' 
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' 
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' 
#' @param quantile A float that corresponds to the quantile used to align the data.
#' 
#' @return A normalized numeric matrix
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[['original_log']]
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- QuantileCentering(assay(obj), conds, type="within conditions", subset.norm=1:10)
#' 
#' @export
#' 
#' @importFrom stats quantile
#' 
QuantileCentering <- function(qData, 
                              conds=NULL,
                              type="overall", 
                              subset.norm=NULL, 
                              quantile=0.15){

  if( missing(conds))
    stop("'conds' is required")
  if( missing(type))
    stop("'type' is required")
  
  
  if (!(type %in% c('overall', 'within conditions')))
    stop("'type' must have one of the following values: 'overall', 'within conditions'")
  
  
  qData <- as.matrix(qData)
  
  if(is.null(subset.norm) || length(subset.norm) < 1)
    subset.norm=1:nrow(qData)

  
  q <- function(x) { 
    stats::quantile(x, probs=quantile, na.rm=TRUE) 
    }
  
  if(length(subset.norm) == 1)
    quantileOverSamples=qData[subset.norm,]
  else
    quantileOverSamples <- apply(qData[subset.norm,], 2, q)

  
  if (type == "overall"){
    cOverall <- q(quantileOverSamples)
    qData <- sweep(qData, 2, quantileOverSamples)
    qData <- qData + cOverall
  } else if (type == "within conditions"){
    qData <- sweep(qData, 2, quantileOverSamples)
    cCond <- NULL
    for (l in unique(conds)) {
      indices <- which(conds== l)
      cCond[l] <- q(quantileOverSamples[indices])
      qData[,indices] <- qData[,indices] + cCond[l]
    }
  }
  return(qData)
}


#' @title Normalisation MeanCentering
#' 
#' @param qData xxx
#' 
#' @param conds xxx
#' 
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' 
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' 
#' @param scaling A boolean that indicates if the variance of the data have to
#' be forced to unit (variance reduction) or not.
#' 
#' @return A normalized numeric matrix
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original_log']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- MeanCentering(qData, conds, type="overall")
#' 
#' @export
#' 
MeanCentering <- function(qData, 
                          conds, 
                          type='overall', 
                          subset.norm=NULL, 
                          scaling=FALSE) {

  if( missing(conds))
    stop("'conds' is required")
  
  qData <- as.matrix(qData)
  
  if(length(subset.norm) == 1)
    meanOverSamples = qData[subset.norm,]
  else
    meanOverSamples <- apply(qData[subset.norm,], 2, mean, na.rm = TRUE)

  
  if (type == "overall"){
    cOverall <- mean(meanOverSamples)
    qData <- sweep(qData, 2, meanOverSamples)
    if (scaling){
      qData <- scale(qData,center=FALSE,scale=TRUE)
      attr(qData,"scaled:scale")<-NULL
    }
    qData <- qData + cOverall
  }
  else if (type == "within conditions"){
    .temp <- sweep(qData, 2, meanOverSamples)
    if (scaling){
      qData <- scale(qData,center=FALSE, scale=TRUE)
      attr(qData,"scaled:scale")<-NULL
    }
    cCond <- NULL
    for (l in unique(conds)) {
      indices <- which(conds== l)
      cCond[l] <- mean(meanOverSamples[indices])
      qData[,indices] <- .temp[,indices] + cCond[l]
    }
  }
  return(qData)
}


#' @title Normalisation vsn
#' 
#' @param qData A numeric matrix.
#' 
#' @param conds xxx
#' 
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' 
#' @return A normalized numeric matrix
#' 
#' @author Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original_log']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- vsn(qData, conds, type="overall")
#' 
#' @export
#' 
#' @importFrom vsn vsnMatrix predict
#' 
vsn = function(qData, 
               conds, 
               type=NULL) {
  if( missing(conds))
    stop("'conds' is required")
  
  if(type == "overall"){
    vsn.fit <- vsn::vsnMatrix(2^qData)
    qData <- vsn::predict(vsn.fit, 2^qData)
  } else if(type == "within conditions"){
    for (l in unique(conds)) {
      indices <- which(conds == l)
      vsn.fit <- vsn::vsnMatrix(2^(qData[,indices]))
      qData[,indices] <- vsn::predict(vsn.fit, 2^qData[,indices])
    }
  }
  return(qData)
}


#' @title Normalisation LOESS
#' 
#' @param qData A numeric matrix.
#' 
#' @param conds xxx
#' 
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' 
#' @param span xxx
#' 
#' @return A normalized numeric matrix
#' 
#' @author Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original_log']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- LOESS(qData, conds, type="overall")
#' 
#' @importFrom limma normalizeCyclicLoess
#' 
#' @export
#' 
LOESS <- function(qData, 
                  conds, 
                  type='overall', 
                  span=0.7) {
  if( missing(conds))
    stop("'conds' is required")
  
  if(type == "overall")
    qData <- limma::normalizeCyclicLoess(x = qData, 
                                         method = "fast", 
                                         span = span)
  else if(type == "within conditions")
    for (l in unique(conds)) {
      indices <- which(conds == l)
      qData[,indices] <- limma::normalizeCyclicLoess(x = qData[,indices],
                                                     method = "fast", 
                                                     span = span)
    }

  return(qData)
}