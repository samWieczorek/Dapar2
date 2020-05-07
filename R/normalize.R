#' Provides several methods to normalize quantitative data from a numeric matrix.
#' They are organized in six main families : GlobalQuantileAlignement, sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
#' For the first family, there is no type.
#' For the five other families, two type categories are available : 
#' "Overall" which means that the value for each protein (ie row) is computed over all the samples ;
#' "within conditions" which means that the value for each protein (ie row) is computed condition by condition.
#' Normalization can be based only on a subset of protein.
#' Proteins with NA in this subset are ignored.
#'
#' @title Normalisation
#' @param obj A \code{SummerizedExperiment} object
#' @param conds xxx
#' @param method "GlobalQuantileAlignment" (for normalizations of important magnitude),
#' "SumByColumns", "QuantileCentering", "Mean Centering", "LOESS" and "vsn".
#' @param type "overall" (shift all the sample distributions at once) or "within conditions" (shift the sample
#' distributions within each condition at a time).
#' @param scaling A boolean that indicates if the variance of the data have to
#' be forced to unit (variance reduction) or not.
#' @param quantile A float that corresponds to the quantile used to align the
#' data.
#' @param span Parameter for LOESS method
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' @return An instance of class \code{SummerizedExperiment} where the quantitative
#' data in the \code{array()} tab has been normalized.
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[['original']]
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' r<-wrapper.normalizeD(obj, "GlobalQuantileAlignment")
#' r<-wrapper.normalizeD(obj, "SumByColumns", conds, type="within conditions")
#' r<-wrapper.normalizeD(obj, "QuantileCentering", conds, type="within conditions")
#' r<-wrapper.normalizeD(obj, "MeanCentering", conds, type="overall")
#' r<-wrapper.normalizeD(obj, "vsn", conds, type="within conditions")
#' r<-wrapper.normalizeD(obj, "LOESS", conds, type="within conditions")
#' @export

wrapper.normalizeD <- function(obj, method, 
                               conds=NULL, type=NULL, scaling=FALSE, quantile=0.15, span = 0.7,
                               subset.norm=NULL){
  
  parammethod<-c("GlobalQuantileAlignment",
                 "SumByColumns",
                 "QuantileCentering",
                 "MeanCentering",
                 "LOESS",
                 "vsn")
  
  if (sum(is.na(match(method, parammethod)==TRUE))>0){
    warning("Parameter method is not correct")
    return (NULL)
  }
  
  paramtype<-c(NULL,"overall", "within conditions")
  if (sum(is.na(match(type, paramtype)==TRUE))>0){
    warning("Parameter type is not correct")
    return (NULL)
  }
  
  metadata(obj)$params$normalizationMethod <- method
  if (method != "GlobalQuantileAlignment") { metadata(obj)$params$normalizationType <- type }
  
  qData <- assay(obj)
  #msg_method <- paste("Normalisation using method =", method,  sep="")
  #msg_type <- paste("With type =", type,  sep="")
  
  if(!is.null(subset.norm)){
    subset.norm=subset.norm[complete.cases(qData)[subset.norm]]
  }
  if(is.null(subset.norm) | length(subset.norm)<1){
    subset.norm=1:nrow(qData)
  }
  
  switch(method,
         GlobalQuantileAlignment = {
           e <- GlobalQuantileAlignment(qData)
           #obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
         },
         SumByColumns = {
           e <- SumByColumns(qData, conds, type, subset.norm)
           #obj@processingData@processing <- c(qData@processingData@processing, msg_method, msg_type)
         },
         QuantileCentering = {
           e <- QuantileCentering(qData, conds, type, subset.norm, quantile)
           metadata(obj)$params$normalizationQuantile <- quantile
           #msg_quantile <- paste("With quantile =", quantile,  sep="")
           #obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_quantile)
         },
         MeanCentering = {
           e <- MeanCentering(qData,  conds, type, subset.norm, scaling)
           metadata(obj)$params$normalizationScaling <- scaling
           #msg_scaling <- paste("With scaling =", scaling,  sep="")
           #obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_scaling)
         },
         vsn = {
           e <- vsn(qData, conds, type)
           #obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
         },
         LOESS = {
           e <- LOESS(qData, conds, type, span)
           metadata(obj)$params$normalizationSpan <- span
           #msg_loess <- paste("With span =", span,  sep="")
           #obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_loess)
         }
  )
  rownames(e) <- rownames(qData)
  colnames(e) <- colnames(qData)
  assay(obj) <- e
  
  return(obj)
}


#' @title Normalisation GlobalQuantileAlignement
#' @param qData A numeric matrix.
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' normalized <- GlobalQuantileAlignment(qData)
#' @export
#' @importFrom preprocessCore normalize.quantiles
GlobalQuantileAlignment <- function(qData) {
  e <- preprocessCore::normalize.quantiles(qData)
  return(e)
}


#' @title Normalisation SumByColumns
#' @param qData A numeric matrix.
#' @param conds xxx
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' conds <- colData(Exp1_R25_pept)[["Condition"]]
#' normalized <- SumByColumns(qData, conds, type="within conditions", subset.norm=1:10)
#' @export
#' @importFrom stats median
SumByColumns <- function(qData, conds=NULL, type=NULL, subset.norm=NULL) {

  e <- 2^qData
  if (type == "overall"){
    if(length(subset.norm)==1){
      sums_cols=e[subset.norm,]
    }else{
      sums_cols <- colSums(e[subset.norm,], na.rm=TRUE)
    }
    for ( i in 1:nrow(e)) {
      e[i, ] <- (e[i, ] / sums_cols)*(stats::median(sums_cols))
    }
  } else if (type == "within conditions"){
    for (l in unique(conds)) {
      indices <- which(conds== l)
      if(length(subset.norm)==1){
        sums_cols=e[subset.norm,indices]
      }else{
        sums_cols <- colSums(e[subset.norm,indices], na.rm=TRUE)
      }
      for (i in 1:nrow(e)){
        e[i,indices] <- (e[i,indices]/sums_cols) * stats::median(sums_cols)
      }
    }
  }
  e <- log2(e)
  return(e)
}


#' @title Normalisation QuantileCentering
#' @param qData A numeric matrix.
#' @param conds xxx
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' @param quantile A float that corresponds to the quantile used to align the data.
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- QuantileCentering(qData, conds, type="within conditions", subset.norm=1:10)
#' @export
#' @importFrom stats quantile
QuantileCentering <- function(qData, conds=NULL, type=NULL, subset.norm=NULL, quantile=0.15){

  q <- function(x) { stats::quantile(x, probs=quantile, na.rm=TRUE) }
  
  if(length(subset.norm)==1){
    quantileOverSamples=qData[subset.norm,]
  }else{
    quantileOverSamples <- apply(qData[subset.norm,], 2, q)
  }
  
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
#' @param qData A numeric matrix.
#' @param conds xxx
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' @param scaling A boolean that indicates if the variance of the data have to
#' be forced to unit (variance reduction) or not.
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- MeanCentering(qData, conds, type="overall")
#' @export
MeanCentering <- function(qData, conds=NULL, type=NULL, subset.norm=NULL, scaling=FALSE) {

  if(length(subset.norm)==1){
    meanOverSamples=qData[subset.norm,]
  } else{
    meanOverSamples <- apply(qData[subset.norm,], 2, mean, na.rm = TRUE)
  }
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
      #attr(qData),"scaled:scale")<-NULL
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
#' @param qData A numeric matrix.
#' @param conds xxx
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- vsn(qData, conds, type="overall")
#' @export
vsn = function(qData, conds, type=NULL) {
  if(type == "overall"){
    vsn.fit <- vsn::vsnMatrix(2^(qData))
    qData <- vsn::predict(vsn.fit, 2^(qData))
  } else if(type == "within conditions"){
    for (l in unique(conds)) {
      indices <- which(conds == l)
      vsn.fit <- vsn::vsnMatrix(2^(qData[,indices]))
      qData[,indices] <- vsn::predict(vsn.fit, 2^(qData[,indices]))
    }
  }
  return(qData)
}


#' @title Normalisation LOESS
#' @param qData A numeric matrix.
#' @param conds xxx
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition at a time).
#' @param span xxx
#' @return A normalized numeric matrix
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[['original']])
#' conds <- colData(Exp1_R25_pept)[['Condition']]
#' normalized <- LOESS(qData, conds, type="overall")
#' @export
LOESS <- function(qData, conds, type=NULL, span=0.7) {
  if(type == "overall"){
    qData <- limma::normalizeCyclicLoess(x = qData, method = "fast", span = span)
  }else if(type == "within conditions"){
    for (l in unique(conds)) {
      indices <- which(conds == l)
      qData[,indices] <- limma::normalizeCyclicLoess(x = qData[,indices],method = "fast", span = span)
    }
  }
  return(qData)
}