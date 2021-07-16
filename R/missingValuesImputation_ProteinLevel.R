
#' @title xxxxxx
#' 
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- impute_dapar(obj, 2,'foo1',  'POV_det_quant', conds=colData(obj)$Condition)
#' 
"impute_dapar"

#' @param object A `SummarizedExperiment` or `QFeatures` object with
#'     missing values to be imputed.
#'     
#' @param i xxxx
#' 
#' @param name xxx
#'     
#' @param method `character(1)` defining the imputation method. See
#'     `imputeMethods()` for available ones. See
#'     [MsCoreUtils::impute_matrix()] for details.
#'     
#' @param ... Additional parameters passed to the inner imputation
#'     function. See [MsCoreUtils::impute_matrix()] for details.
#'     
#' @export
#' 
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname impute_dapar
#' 
setMethod("impute_dapar", "SummarizedExperiment",
          function(object, method, ...) {
              res <- impute_matrix_dapar(assay(object), method, ...)
              rownames(res) <- rownames(assay(object))
              colnames(res) <- colnames(assay(object))
              SummarizedExperiment::assay(object) <- res
              object
          })


#' @param i Defines which element of the `QFeatures` instance to
#'     impute. If missing, all assays will be imputed.
#' 
#' @export
#' 
#' @importFrom SummarizedExperiment assay
#' 
#' @rdname impute_dapar
#' 
setMethod("impute_dapar", "QFeatures",
          function(object, i, name, method, ...) {
              if (missing(i))
                  i  <-  length(object)
              if (is.numeric(i)) i <- names(object)[[i]]
              
              argg <- c(as.list(environment()), list(...))
              
              tmp <- object[[i]]
              res <- impute_matrix_dapar(assay(tmp), method, ...)
              SummarizedExperiment::assay(tmp) <- res
              #metadata(tmp)$Params <-  argg[-match(c('object', 'i', 'name'), names(argg))]
              p <- argg[-match(c('object', 'i', 'name'), names(argg))]
              metadata(tmp)$Params <-  paste0(names(p),'=', p, collapse=', ')
              object <- QFeatures::addAssay(object,
                                           tmp,
                                           name)
              QFeatures::addAssayLinkOneToOne(object, from = i, to = name)
          })




#' @title xxx
#' 
#' @param x xxxx
#' 
#' @param method xxxx
#' 
#' @param ... xxxxx
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' imp <- impute_matrix_dapar(assay(obj[[2]]), method='knn_by_conds', colData(obj)$Condition, 3)
#' 
#' imp <- impute_matrix_dapar(assay(obj[[2]]), method='pa', colData(obj)$Condition, q.min=0.03)
#' 
#' imp <- impute_matrix_dapar(assay(obj[[2]]), method='det_quant')
#' 
#' imp <- impute_matrix_dapar(assay(obj[[2]]), method='slsa', sampleTab = colData(obj))
#' 
#' @export
#' 
#' @rdname impute_matrix_dapar
#' 
impute_matrix_dapar <- function(x,
                          method,
                          ...) {
    if (!anyNA(x)) return(x)
    if (missing(method))
        stop("Please specify an imputation method. ")
  
    if (!(method %in% imputeMethodsDapar()))
        stop("'method' is not a valid method in DAPAR.")
    method <- match.arg(method,
                        choices = imputeMethodsDapar(),
                        several.ok = FALSE)
    res <- x
    
    if (method == "POV_knn_by_conds") {
      res <- POV_impute_knn_by_conditions(x, ...)
    } else if (method == "knn_by_conds") {
        res <- impute_knn_by_conditions(x, ...)
    } else if (method == "pa") {
        requireNamespace("imp4p")
        res <- impute_pa(x, ...)
    } else if (method == "det_quant") {
        res <- impute_det_quant(x, ...)
    } else if (method == "POV_det_quant") {
      res <- POV_impute_det_quant(x, ...)
    }else if (method == "slsa") {
        res <- impute_slsa(x, ...)
    } else if (method == "POV_slsa") {
      res <- POV_impute_slsa(x, ...)
    } else if (method == "mle_dapar") {
        res <- impute_mle_dapar(x, ...)
    } else if (method == "mi") {
        res <- impute_mi(x, ...)
    } else if (method == "fixed_val") {
      res <- impute_fixed_value(x, ...)
    }
    ## else method == "none" -- do nothing
    res
}



#' @title xxxx
#' 
#' @export
#' 
imputeMethodsDapar <- function()
    c('POV_knn_by_conds', "knn_by_conds", "pa", "det_quant", "slsa", "mle_dapar", 'fixed_val', 'mi', 'POV_slsa', 'POV_det_quant', "none")

#' @title List the methods available in DAPAR to impute Partially Observed Values (POV)
#' 
#' @export
#' 
impute_POV_Methods <- function()
  c(imputeMethodsDapar()[grep('POV', imputeMethodsDapar())],"none")


#' @title List the methods in DAPAR to impute Missing Entire Condition (MEC)
#' 
#' @export
#' 
impute_MEC_Methods <- function()
  c("knn_by_conds", "pa", "det_quant", "slsa", "mle_dapar", "none")


#' @title KNN missing values imputation by conditions
#' 
#' @description This function imputes the missing values condition by condition.
#' 
#' @param x xxxx
#' 
#' @param conds xxx
#' 
#' @param k the number of neighbors.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' imp <- POV_impute_knn_by_conditions(assay(obj[[2]]), colData(obj)$Condition, 3)
#' 
#' @export
#' 
#' @importFrom impute impute.knn
#' 
POV_impute_knn_by_conditions <- function(x, conds=NULL, k=3){
  
  if(missing(conds))
    stop("'conds' is missing.")
  
  if (is.null(conds))
    stop("'conds' is empty.")
  
  MECIndex <- find_MEC_matrix(x, conds)

  res <-x
  u_conds <- unique(conds)
  
  
  for (i in 1:length(u_conds)){
    ind <- which(conds == u_conds[i])
    resKNN <- impute::impute.knn(res[,ind] ,k = k, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000,1))
    res[,ind] <- resKNN[[1]]
  }
  
  res <- restore_MEC_matrix(res, conds, MECIndex)
  
  return(res)
}



#' @title KNN missing values imputation by conditions
#' 
#' @description This function imputes the missing values condition by condition.
#' 
#' @param x xxxx
#' 
#' @param conds xxx
#' 
#' @param k the number of neighbors.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' imp <- impute_knn_by_conditions(assay(obj[[2]]), colData(obj)$Condition, 3)
#' 
#' @export
#' 
#' @importFrom impute impute.knn
#' 
impute_knn_by_conditions <- function(x, conds=NULL, k=3){
    
    if(missing(conds))
       stop("'conds' is missing.")

    if (is.null(conds))
        stop("'conds' is empty.")

    res <-x
    u_conds <- unique(conds)
    
    
    for (i in 1:length(u_conds)){
        ind <- which(conds == u_conds[i])
        resKNN <- impute::impute.knn(res[,ind] ,k = k, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000,1))
        res[,ind] <- resKNN[[1]]
    }

    return(res)
}



#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param x xxxxx
#' 
#' @param conds xxxx
#' 
#' @param q.min Same as the function \code{impute.pa} in the package \code{imp4p}
#' 
#' @return The xxxx.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0.2), name='filtered')
#' imp <- impute_pa(assay(object[['filtered']]), colData(object)$Condition)
#' 
#' @export
#' 
#' @importFrom imp4p impute.pa
#' 
impute_pa <- function(x, conds, q.min = 0.025){
    res <- x
    res <- imp4p::impute.pa(x, 
              conditions=as.factor(conds), 
              q.min = q.min, 
              q.norm = 3,  
              eps=0
              )
    
    return (res[["tab.imp"]])
}





#' @title xxx
#' 
#' @param x An instance of class \code{DataFrame}
#' 
#' @param ... Additional arguments
#' 
#' @return An imputed instance of class \code{DataFrame}
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' impute_det_quant(assay(Exp1_R25_pept[[2]]))
#' 
#' @export
impute_det_quant <- function(x, ...){
    if (is.null(x)){return(NULL)}
    
    values <- getQuantile4Imp(x, ...)
    res <- x
   for(i in 1:dim(x)[2]){
        col <- x[,i]
        col[which(is.na(col))] <- values$shiftedImpVal[i]
        x[,i] <- col
   }
   
    return(x)
}

#' @title Imputation of the POV with the method of deterministic Quantile
#' 
#' @param x An instance of class \code{DataFrame}
#' 
#' @param conds xxx
#' 
#' @param ... Additional arguments
#' 
#' @return An imputed instance of class \code{DataFrame}
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' POV_impute_det_quant(assay(Exp1_R25_pept[[2]]), colData(Exp1_R25_pept)$Condition)
#' 
#' @export
POV_impute_det_quant <- function(x, conds, ...){
  if (is.null(x)){return(NULL)}
  
  MECIndex <- find_MEC_matrix(x, conds)
  
  values <- getQuantile4Imp(x, ...)
  res <- x
  for(i in 1:dim(x)[2]){
    col <- x[,i]
    col[which(is.na(col))] <- values$shiftedImpVal[i]
    x[,i] <- col
  }
  
  x <- restore_MEC_matrix(x, conds, MECIndex)
  
  return(x)
}

#' @title Quantile imputation value definition
#' 
#' @description This method returns the q-th quantile of each column of an expression set, up to a scaling factor
#' 
#' @param x An expression set containing quantitative values of various replicates
#' 
#' @param qval The quantile used to define the imputation value
#' 
#' @param factor A scaling factor to multiply the imputation value with
#' 
#' @return A list of two vectors, respectively containing the imputation values and the rescaled imputation values
#' 
#' @author Thomas Burger
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' qData <- assay(Exp1_R25_pept[[2]])
#' getQuantile4Imp(qData) 
#' 
#' @export
#' 
#' @importFrom stats quantile
#' 
getQuantile4Imp <- function(x, qval=0.025, factor=1){
    r1 <- apply(x, 2, quantile, qval, na.rm=TRUE)
    r2 <- r1*factor
    return(list(ImpVal = r1, shiftedImpVal = r2))
}




#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param x A DataFrame.
#' 
#' @param sampleTab xxxxx
#' 
#' @return A DataFrame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' imp <- impute_slsa(assay(obj[[2]]), colData(obj))
#' 
#' @export
#' 
#' @importFrom imp4p impute.slsa
#' 
impute_slsa <- function(x, sampleTab){
    # sort conditions to be compliant with impute.slsa (from imp4p version 0.9) which only manage ordered samples 
    old.sample.name <- sampleTab$Sample.name
    new.order <- order(sampleTab$Condition)
     new.sampleTab <- sampleTab[new.order,]
    conds <- factor(new.sampleTab$Condition, levels=unique(new.sampleTab$Condition))
    new.x <- x[,new.order]
    
    res <- imp4p::impute.slsa(new.x, conditions=conds, nknn=6, selec="all", weight=1, ind.comp=1)
    #restore old order
    res <- res[,old.sample.name]
    
    return (res)
}



#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param x A DataFrame.
#' 
#' @param sampleTab xxxxx
#' 
#' @return A DataFrame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' imp <- POV_impute_slsa(assay(obj[[2]]), colData(obj))
#' 
#' @export
#' 
#' @importFrom imp4p impute.slsa
#' 
POV_impute_slsa <- function(x, sampleTab){
  
  MECIndex <- find_MEC_matrix(x, sampleTab$Condition)
  
  # sort conditions to be compliant with impute.slsa (from imp4p version 0.9) which only manage ordered samples 
  old.sample.name <- sampleTab$Sample.name
  new.order <- order(sampleTab$Condition)
  new.sampleTab <- sampleTab[new.order,]
  conds <- factor(new.sampleTab$Condition, levels=unique(new.sampleTab$Condition))
  new.x <- x[,new.order]
  
  res <- imp4p::impute.slsa(new.x, conditions=conds, nknn=6, selec="all", weight=1, ind.comp=1)
  #restore old order
  res <- res[,old.sample.name]
  
  res <- restore_MEC_matrix(res, sampleTab$Condition, MECIndex)
  
  return (res)
}


#' This method is a wrapper to
#' objects of class \code{MSnSet} and imputes missing values with a fixed value.
#'
#' @title Missing values imputation from a \code{MSnSet} object
#' 
#' @param x An assay of a SummarizedExperiment
#' 
#' @param value A float .
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' impute_fixed_value(assay(Exp1_R25_pept[1:1000], 2), 0.001)
#' 
#' @export
#' 
impute_fixed_value <- function(x, value){
  
  x[is.na(x)] <- value
  return (x)
}



#' 
#'
#' @title Put back LAPALA into xxxx
#' 
#' @description This method is used to put back the LAPALA that have been identified previously
#' 
#' @param x xxxxx
#' 
#' @param conds xxx
#' 
#' @param MECIndex A data.frame that contains index of MEC (see findMECBlock) .
#' 
#' @return The object \code{obj} where LAPALA have been reintroduced
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' lapala <- find_MEC_matrix(assay(obj[[2]]), colData(obj)$Condition)
#' assay(obj[[2]]) <- impute_det_quant(assay(obj[[2]]))
#' assay(obj[[2]]) <- restore_MEC_matrix(assay(obj[[2]]), colData(obj)$Condition, lapala)
#' 
#' @export
#' 
restore_MEC_matrix <- function(x, conds, MECIndex){
    u_conds <- unique(conds)
    for (i in 1:nrow(MECIndex))
    {
       replicates <- which(conds == u_conds[MECIndex[i,"Condition"]])
        x[MECIndex[i,"Line"], as.vector(replicates)] <- NA
    }
    return(x)
}



#' @title Finds the LAPALA into a \code{MSnSet} object
#' 
#' @param x xxxxx
#' 
#' @param conds xxx
#' 
#' @return A data.frame that contains the indexes of LAPALA
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' lapala <- find_MEC_matrix(assay(obj[[2]]), colData(obj)$Condition)
#' 
#' @export
#' 
find_MEC_matrix <- function(x, conds){
    
 # if(class(x)!='matrix')
 #   stop("'c' must be an instance of class 'SummarizedExperiment'.")
  
  
    u_conds <- unique(conds)
    nbCond <- length(u_conds)
    
    s <- data.frame()
    
    for (i in 1:length(u_conds)){
        ind <- which(conds == u_conds[i])
        lNA <- which(apply(is.na(x[,ind]), 1, sum)==length(ind))
        if (length(lNA) > 0)
        {
            tmp <- data.frame(i,which(apply(is.na(x[,ind]), 1, sum)==length(ind)))
            names(tmp) <- c("Condition", "Line")
            s <- rbind(s,tmp)
        }
    }
    return(s)
}



























