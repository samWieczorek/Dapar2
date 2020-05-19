#' The function impute_fixed_value is deleted and the dev must use the function impite_with from the package MSCoreUtils



#' @param object A `SummarizedExperiment` or `Features` object with
#'     missing values to be imputed.
#' @param method `character(1)` defining the imputation method. See
#'     `imputeMethods()` for available ones. See
#'     [MsCoreUtils::impute_matrix()] for details.
#' @param ... Additional parameters passed to the inner imputation
#'     function. See [MsCoreUtils::impute_matrix()] for details.
#' 
#' @export
#' @rdname impute
setMethod("impute_dapar", "SummarizedExperiment",
          function(object, method, ...) {
              res <- impute_matrix_dapar(assay(object), method, ...)
              assay(object) <- res
              object
          })


#' @param i Defines which element of the `Features` instance to
#'     impute. If missing, all assays will be imputed.
#' 
#' @export
#' @rdname impute
setMethod("impute_dapar", "Features",
          function(object, method, ..., i) {
              if (missing(i))
                  i  <-  seq_len(length(object))
              for (ii in i) {
                  res <- impute_matrix_dapar(assay(object[[ii]]), method, ...)
                  assay(object[[ii]]) <- res
              }
              object
          })




#' @examples
#' library(Features)
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
impute_matrix_dapar <- function(x,
                          method,
                          ...) {
    if (!anyNA(x)) return(x)
    if (missing(method))
        stop("Please specify an imputation method. ",
             "See '?impute_matrix_dapar' for details.")
    method <- match.arg(method,
                        choices = imputeMethodsDapar(),
                        several.ok = FALSE)
    res <- x
    if (method %in% c("impute_pa"))
        requireNamespace("imp4p")
    
    if (method == "knn_by_conds") {
        res <- impute_knn_by_conditions(x, ...)
    } else if (method == "pa") {
        res <- impute_pa(x, ...)
    } else if (method == "det_quant") {
        res <- impute_det_quant(x, ...)
    } else if (method == "slsa") {
        res <- impute_slsa(x, ...)
    }
    ## else method == "none" -- do nothing
    res
}



##' @export
##' @rdname imputation
imputeMethodsDapar <- function()
    c("knn_by_conds", "pa", "det_quant", "slsa", "none")



#' @title KNN missing values imputation by conditions
#' 
#' @description This function imputes the missing values condition by condition.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param k the number of neighbors.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
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
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param q.min Same as the function \code{impute.pa} in the package \code{imp4p}
#' 
#' @return The \code{exprs(obj)} matrix with imputed values instead of missing values.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(Features)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, Features::filterNA(object[[2]],  pNA = 0.2), name='filtered')
#' imp <- impute_pa(assay(object[['filtered']]), colData(object)$Condition)
#' 
#' @export
#' 
#' importFrom imp4p impute.pa
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





#' @title Wrapper of the function \code{\link{impute.detQuant}} for objects
#' of class \code{MSnSet}
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
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' impute_det_quant(assay(Exp1_R25_pept[[2]])
#' 
#' @export
impute_det_quant <- function(x,...){
    if (is.null(obj)){return(NULL)}
    
    values <- getQuantile4Imp(x, ...)
    res <- x
   for(i in 1:dim(x)[2]){
        col <- x[,i]
        col[which(is.na(col))] <- values$shiftedImpVal[i]
        x[,i] <- col
    }
    return(x)
}



#' @title Quantile imputation value definition
#' 
#' @description This method returns the q-th quantile of each colum of an expression set, up to a scaling factor
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
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
#' @param conds xxxxx
#' 
#' @return A DataFrame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
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




#' 
#'
#' @title Put back LAPALA into  a \code{MSnSet} object
#' 
#' @description This method is used to put back the LAPALA that have been identified previously
#' 
#' @param obj An object of class \code{Su}.
#' 
#' @param MECIndex A data.frame that contains index of MEC (see findMECBlock) .
#' 
#' @return The object \code{obj} where LAPALA have been reintroduced
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' lapala <- find_MEC(assay(obj[[2]]), colData(obj)$Condition)
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
#' @param obj An object of class \code{SummarizedExperiment}.
#' 
#' @return A data.frame that contains the indexes of LAPALA
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' lapala <- find_MEC_matrix(assay(obj[[2]]), colData(obj)$Condition)
#' 
#' @export
#' 
find_MEC_matrix <- function(x, conds){
    
    u_conds <- unique(conds)
    nbCond <- length(conditions)
    
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



























