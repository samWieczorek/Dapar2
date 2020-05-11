

setGeneric("computeTTest", function(object, ...) standardGeneric("computeTTest"))

##' @title Features compute t-tests
##'
##' @description
##'
##' This manual page describes the computation of statistical test using [Features] objects. In the following
##' functions, if `object` is of class `Features`, and optional assay
##' index or name `i` can be specified to define the assay (by name of
##' index) on which to operate.
##'
##' The following functions are currently available:
##'
##' - `compute.t.test(object, base = 2, i, pc = 0)` log-transforms (with
##'   an optional pseudocount offset) the assay(s).
##'
##' - `compute.group.t.test(object, method, i)` normalises the assay(s) according
##'   to `method` (see Details).
##' 
##'
##' @details
##' 
##' The `method` parameter in `normalize` can be one of `"sum"`,
##' `"max"`, `"center.mean"`, `"center.median"`, `"div.mean"`,
##' `"div.median"`, `"diff.meda"`, `"quantiles`", `"quantiles.robust`"
##' or `"vsn"`. The [MsCoreUtils::normalizeMethods()] function
##' returns a vector of available normalisation methods.
##'
##' - For `"sum"` and `"max"`, each feature's intensity is divided by the
##'   maximum or the sum of the feature respectively. These two methods are
##'   applied along the features (rows).
##'
##' - `"center.mean"` and `"center.median"` center the respective sample
##'   (column) intensities by subtracting the respective column means or
##'   medians. `"div.mean"` and `"div.median"` divide by the column means or
##'   medians.
##'
##' - `"diff.median"` centers all samples (columns) so that they all match the
##'   grand median by subtracting the respective columns medians differences to
##'   the grand median.
##'
##' - Using `"quantiles"` or `"quantiles.robust"` applies (robust) quantile
##'   normalisation, as implemented in [preprocessCore::normalize.quantiles()]
##'   and [preprocessCore::normalize.quantiles.robust()]. `"vsn"` uses the
##'   [vsn::vsn2()] function.  Note that the latter also glog-transforms the
##'   intensities.  See respective manuals for more details and function
##'   arguments.
##'
##' For further details and examples about normalisation, see
##' [MsCoreUtils::normalize_matrix()].
##'
##' @param  object An object of class `Features` or `SummarizedExperiment`.
##'
##' @param base `numeric(1)` providing the base with respect to which
##'     logarithms are computed. Defaults is 2.
##'
##' @param pc `numeric(1)` with a pseudocount to add to the
##'     quantitative data. Useful when (true) 0 are present in the
##'     data. Default is 0 (no effect).
##'
##' @param center `logical(1)` (default is `TRUE`) value or
##'     numeric-alike vector of length equal to the number of columns
##'     of `object`. See [base::scale()] for details.
##' 
##' @param scale `logical(1)` (default is `TRUE`) or a numeric-alike
##'     vector of length equal to the number of columns of
##'     `object`. See [base::scale()] for details.
##'
##' @param method `character(1)` defining the normalisation method to
##'     apply. See Details.
##' 
##' @param i A numeric vector or a character vector giving the index or the 
##'     name, respectively, of the assay(s) to be processed.
##'
##' @param name A `character(1)` naming the new assay name. Defaults
##'     are `logAssay` for `logTransform`, `scaledAssay` for
##'     `scaleTranform` and `normAssay` for `normalize`.
##'
##' @param ... Additional parameters passed to inner functions.
##'
##' @aliases logTransform logTransform,SummarizedExperiment-method logTransform,Features-method
##'
##' @aliases scaleTransform scaleTransform,SummarizedExperiment-method scaleTransform,Features-method
##' 
##' @aliases normalize normalize,SummarizedExperiment-method normalize,Features-method
##'
##' @aliases normalizeMethods
##'
##' @name Features-processing
##'
##' @rdname Features-processing
##'
##' @examples
##'
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept
#' object <- addAssay(object, Features::filterNA(object[[2]],  pNA = 0), name='filtered')
#' sTab <- colData(object)
#' ttest <- computeTTest(object[[3]], sTab, contrast="OnevsOne")
#' 
#' object <-computeTTest(object, 3, name = "ttestAssay", contrast = 'OnevsOne')
#' 
NULL

##' @exportMethod computeTTest
 #' @rdname Features-dapar-compute-ttests
setMethod("computeTTest", "SummarizedExperiment",
          function(object,
                   sampleTab,
                   ...) {
              df <- compute.t.test(assay(object), sampleTab, ...)
              #rownames(df) <- rownames(assay(object))
              e <- object
              metadata(e)$t_test <- df
              e
          })


#' @rdname Features-dapar-compute-ttests
setMethod("computeTTest", "Features",
          function(object, i, name = "ttestAssay", ...) {
              if (missing(i))
                  stop("Provide index or name of assay to be processed")
              if (length(i) != 1)
                  stop("Only one assay to be processed at a time")
              if (is.numeric(i)) i <- names(object)[[i]]
              object <- addAssay(object,
                                 computeTTest(object[[i]], sampleTab = colData(object), ...),
                                 name)
              addAssayLinkOneToOne(object, from = i, to = name)
          })




#' This function is xxxxxx
#'
#' @title xxxxxx
#' 
#' @description 
#' 
#' @param qData A matrix of quantitative data, without any missing values.
#' 
#' @param sampleTab A data.frame which contains the samples data.  
#' 
#' @param contrast Indicates if the test consists of the comparison of each 
#' biological condition versus 
#' each of the other ones (contrast=1; 
#' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#' 
#' @param type A string which is the type of statistical test to proceed. Available values are 'Student' and 'Welch'.
#' Default value is 'Student'.
#' 
#' @return A list of two items : logFC and P_Value; both are dataframe. The first one contains
#' the logFC values of all the comparisons (one column for one comparison), the second one contains
#' the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes
#' are identical and correspond to the description of the comparison. 
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' Exp1_R25_pept <- addAssay(Exp1_R25_pept, Features::filterNA(Exp1_R25_pept[[2]],  pNA = 0), name='filtered')
#' sTab <- as.data.frame(colData(Exp1_R25_pept))
#' qData <- assay(Exp1_R25_pept[['filtered']])
#' ttest <- compute.t.test(qData, sTab ,"OnevsOne")
#' 
#' @export
#' 
#' @importFrom stats t.test
#' 
compute.t.test <- function(qData, sampleTab, contrast="OnevsOne", type="Student"){

    
    # if (class(sampleTab) != 'data.frame'){
    #     stop("'sampleTab' is not of class data.frame.")
    # }
    
    .type <- type=='Student'
    sampleTab <- as.data.frame(sampleTab)
    
res<-list()
logFC <- list()
P_Value <- list()

nbComp <- NULL

sampleTab.old <- sampleTab
Conditions.f <- factor(sampleTab$Condition, levels=unique(sampleTab$Condition))
sampleTab <- sampleTab[unlist(lapply(split(sampleTab, Conditions.f), function(x) {x['Sample.name']})),]
qData <- qData[,unlist(lapply(split(sampleTab.old, Conditions.f), function(x) {x['Sample.name']}))]
Conditions <- sampleTab$Condition

Cond.Nb<-length(levels(Conditions.f))


    if(contrast=="OnevsOne"){
        nbComp <- Cond.Nb*(Cond.Nb-1)/2

        for(i in 1:(Cond.Nb-1)){
            for (j in (i+1):Cond.Nb){
    
                c1Indice <- which(Conditions==levels(Conditions.f)[i])
                c2Indice <- which(Conditions==levels(Conditions.f)[j])
    
                res.tmp <- apply(qData[,c(c1Indice,c2Indice)], 1, 
                                 function(x) {
                   stats::t.test(x~Conditions[c(c1Indice,c2Indice)],  var.equal=.type)
                })
                p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
                m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
                m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
                m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
                m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
                logFC.tmp <- m1.tmp - m2.tmp
                if (grepl(levels(Conditions.f)[i], m2.name)){logFC.tmp <- -logFC.tmp}
                
                txt <- paste(levels(Conditions.f)[i],"_vs_",levels(Conditions.f)[j], sep="")
                
                logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
                P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
            }
        }
    } ##end Contrast==1

    if(contrast=="OnevsAll"){
        nbComp <- Cond.Nb
        
        for(i in 1:nbComp){
            
            c1 <- which(Conditions==levels(Conditions.f)[i])
           
            Cond.t.all<-c(1:length(Conditions))
            Cond.t.all[c1]<-levels(Conditions.f)[i]
            Cond.t.all[-c1]<-"all"
            
            res.tmp <- apply(qData, 1, 
                             function(x) {
                                 stats::t.test(x~Cond.t.all, var.equal=.type)
                             })
            
            p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
            m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
            m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
            m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
            m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
            logFC.tmp <- m1.tmp - m2.tmp
            if (grepl(levels(Conditions.f)[i], m2.name)){logFC.tmp <- -logFC.tmp}
            
            txt <- paste(levels(Conditions.f)[i],"_vs_(all-",levels(Conditions.f)[i],")", sep="")
            logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
            P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
        }
    } # End Contrast=2
    

    res.l <- DataFrame(logFC, P_Value)
    colnames(res.l) <- c(names(logFC), names(P_Value))
    
    return(res.l) 
    
}







#' #' This function is xxxxxx
#' #'
#' #' @title xxxxxx
#' #' @param qData A matrix of quantitative data, without any missing values.
#' #' @param sTab xxxx  
#' #' @param contrast Indicates if the test consists of the comparison of each 
#' #' biological condition versus 
#' #' each of the other ones (contrast=1; 
#' #' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#' #' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#' #' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#' #' @param type xxxxx
#' #' @return A list of two items : logFC and P_Value; both are dataframe. The first one contains
#' #' the logFC values of all the comparisons (one column for one comparison), the second one contains
#' #' the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes
#' #' are identical and correspond to the description of the comparison. 
#' #' @author Florence Combes, Samuel Wieczorek
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' #' Exp1_R25_pept <- addAssay(Exp1_R25_pept, Features::filterNA(Exp1_R25_pept[['original_log']],  pNA = 0), name='filtered')
#' #' sTab <- as.data.frame(colData(Exp1_R25_pept)@listData)
#' #' qData <- assay(Exp1_R25_pept[['filtered']])
#' #' ttest <- compute.t.tests2(qData,sTab ,"OnevsOne")
#' #' @export
#' #' @importFrom utils combn
#' #' @importFrom stats t.test
#' compute.t.tests2 <- function(qData,sTab, contrast="OnevsOne", type="Student"){
#'     
#'     
#'     switch(type,
#'            Student=.type <- TRUE,
#'            Welch=.type <- FALSE)
#'     
#'     
#'     
#'     res<-list()
#'     logFC <- list()
#'     P_Value <- list()
#'     
#'     nbComp <- NULL
#'     
#'     sTab.old <- sTab
#'     Conditions.f <- factor(sTab$Condition, levels=unique(sTab$Condition))
#'     sTab <- sTab[unlist(lapply(split(sTab, Conditions.f), function(x) {x['Sample.name']})),]
#'     qData <- qData[,unlist(lapply(split(sTab.old, Conditions.f), function(x) {x['Sample.name']}))]
#'     Conditions <- sTab$Condition
#'     
#'     #Cond.Nb<-length(levels(Conditions.f))
#'     
#'     
#'     if(contrast=="OnevsOne"){
#'         comb <- utils::combn(levels(Conditions.f), 2)
#'         #nbComp <- Cond.Nb*(Cond.Nb-1)/2
#'         
#'         for(i in 1:ncol(comb)){
#'                 
#'                 c1Indice <- which(Conditions==comb[1,i])
#'                 c2Indice <- which(Conditions==comb[2,i])
#'                 
#'                 res.tmp <- apply(qData[,c(c1Indice,c2Indice)], 1, 
#'                                  function(x) {
#'                                      stats::t.test(x~Conditions[c(c1Indice,c2Indice)],  var.equal=.type)
#'                                  })
#'                 p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
#'                 m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
#'                 m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
#'                 m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
#'                 m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
#'                 logFC.tmp <- m1.tmp - m2.tmp
#'                 if (grepl(comb[1,i], m2.name)){logFC.tmp <- -logFC.tmp}
#'                 
#'                 txt <- paste(comb[1,i],"_vs_",comb[2,i], sep="")
#'                 
#'                 logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
#'                 P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
#'             }
#'     } ##end Contrast==1
#'     
#'     if(contrast=="OnevsAll"){
#'         nbComp <- Cond.Nb
#'         
#'         for(i in 1:nbComp){
#'             
#'             c1 <- which(Conditions==levels(Conditions.f)[i])
#'             
#'             Cond.t.all<-c(1:length(Conditions))
#'             Cond.t.all[c1]<-levels(Conditions.f)[i]
#'             Cond.t.all[-c1]<-"all"
#'             
#'             res.tmp <- apply(qData, 1, 
#'                              function(x) {
#'                                  stats::t.test(x~Cond.t.all, var.equal=.type)
#'                              })
#'             
#'             p.tmp <- unlist(lapply(res.tmp,function(x)x$p.value))
#'             m1.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[1])))
#'             m2.tmp <- unlist(lapply(res.tmp,function(x)as.numeric(x$estimate[2])))
#'             m1.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[1])))[1]
#'             m2.name <- names(unlist(lapply(res.tmp,function(x)x$estimate[2])))[1]
#'             logFC.tmp <- m1.tmp - m2.tmp
#'             if (grepl(levels(Conditions.f)[i], m2.name)){logFC.tmp <- -logFC.tmp}
#'             
#'             txt <- paste(levels(Conditions.f)[i],"_vs_(all-",levels(Conditions.f)[i],")", sep="")
#'             logFC[[paste(txt, "logFC", sep="_")]] <- logFC.tmp
#'             P_Value[[paste(txt, "pval", sep="_")]] <- p.tmp
#'         }
#'     } # End Contrast=2
#'     
#'     
#'     res.l <- list(
#'         logFC = as.data.frame(logFC),
#'         P_Value = as.data.frame(P_Value)
#'     )
#'     colnames(res.l$logFC) <- names(logFC)
#'     colnames(res.l$P_Value) <- names(P_Value)
#'     
#'     return(res.l) 
#'     
#' }