
#' @title Check the validity of the experimental design
#'
#' @description
#'
#' This manual page describes the computation of statistical test using [QFeatures] objects. In the following
#' functions, if `object` is of class `QFeatures`, and optional assay
#' index or name `i` can be specified to define the assay (by name of
#' index) on which to operate.
#'
#' The following functions are currently available:
#'
#' - `compute.t.test(xxxxx)` xxxxx.
#'
#' - `compute.group.t.test(xxxxx)` xxxxx.
#'   
#' - `limma.complete.test(object, sampleTab)` uses the package Limma 
#' 
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' sTab <- colData(object)
#' gttest.se <- t_test_sam(object[[3]], sTab, FUN = compute.t.test)
#' object <- addAssay(object, gttest.se, name='t-test')
#' comp <- '25fmol_vs_10fmol'
#' da.se <- diff_analysis_sam(object[['t-test']], comp)
#' 
#' object <- diff_analysis_sam(object, 't-test', name='diffAna', comp)
#' 
"diff_analysis_sam"

#' @param  object An object of class `QFeatures` or `SummarizedExperiment`.
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @rdname diff_analysis_sam
#' 
setMethod("diff_analysis_sam", "SummarizedExperiment",
          function(object,
                   ...) {
            
            res <- diffAnalysis(object, ...)
            res
          })


#' @param  object An object of class `QFeatures` or `SummarizedExperiment`.
#' 
#' @param i A numeric vector or a character vector giving the index or the 
#'     name, respectively, of the assay(s) to be processed.
#'
#' @param name A `character(1)` naming the new assay name. Default is `ttestAssay`.
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @rdname diff_analysis_sam
#' 
setMethod("diff_analysis_sam", "QFeatures",
          function(object, i, name = "diffAnaAssay",  ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            
            
            object <- addAssay(object,
                               diff_analysis_sam(object[[i]],  ...),
                               name)
            addAssayLinkOneToOne(object, from = i, to = name)
          })











#' @title Plots a histogram of p-values 
#' 
#' @param pval_ll A data frame which contains the p-values for the set of comparisons to show.
#' 
#' @param bins xxx
#' 
#' @param pi0 xxx
#' 
#' @return A plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' object <- addListAdjacencyMatrices(object, 3)
#' sTab <- colData(object)
#' se <- t_test_sam(object[[3]], sTab, FUN = compute.t.test)
#' pval <- metadata(se)$t_test[['25fmol_vs_10fmol_pval']]
#' histPValue_HC(pval)
#' 
#' @export
#' 
#' @import highcharter
#' @importFrom graphics hist
#' 
histPValue_HC <- function(pval_ll, bins=80, pi0=1){
  
  h <- graphics::hist(sort(unlist(pval_ll)), freq=F,breaks=bins)
  
  serieInf <- sapply(h$density, function(x)min(pi0, x) )
  serieSup <- sapply(h$density, function(x)max(0, x-pi0) )
  
  hc <- highchart() %>% 
    hc_chart(type = "column") %>%
    hc_add_series(data = serieSup, name="p-value density") %>%
    hc_add_series(data = serieInf, name="p-value density") %>%
    hc_title(text = "P-value histogram") %>% 
    hc_legend(enabled = FALSE) %>%
    hc_colors(c("green", "red")) %>%
    hc_xAxis(title = list(text = "P-value"), categories=h$breaks)%>%
    hc_yAxis(title = list(text="Density"),
             plotLines=list(list(color= "blue" , width = 2, value = pi0, zIndex = 5))) %>%
    hc_tooltip(headerFormat= '',
               pointFormat = "<b> {series.name} </b>: {point.y} ",
               valueDecimals = 2) %>%
    dapar_hc_ExportMenu(filename = "histPVal") %>%
    hc_plotOptions(
      column=list(
        groupPadding= 0,
        pointPadding= 0,
        borderWidth= 0
      ),
      series=list(
        stacking = "normal",
        animation=list(duration = 100),
        connectNulls= TRUE,
        marker=list(enabled = FALSE)
      )
    ) %>%
    hc_add_annotation(
      labelOptions = list(
        backgroundColor = 'transparent',
        verticalAlign = 'top',
        y = -30,
        borderWidth = 0,
        x = 20,
        style=list(
          fontSize= '1.5em',
          color= 'blue'
        )
        
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = 80,
            y = pi0
          ),
          text = paste0("pi0=", pi0)
        )
      )
    )
  return(hc)
}




#' @title Computes the FDR corresponding to the p-values of the differential analysis using
#' 
#' @description This function is a wrappper to the function adjust.p from the cp4p package. It returns the FDR corresponding to the p-values of the 
#' differential analysis. The FDR is computed with the function \code{p.adjust}\{stats\}..
#' 
#' @param logFC The result (logFC values) of the differential analysis processed 
#' by 'limmaCompleteTest'
#' 
#' @param pval The result (p-values) of the differential analysis processed 
#' by 'limmaCompleteTest'
#' 
#' @param threshold_PVal The threshold on p-pvalue to
#' distinguish between differential and non-differential data 
#' 
#' @param threshold_LogFC The threshold on log(Fold Change) to
#' distinguish between differential and non-differential data 
#' 
#' @param pi0Method The parameter pi0.method of the method adjust.p 
#' in the package \code{cp4p}
#' 
#' @return The computed FDR value (floating number)
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' obj <- addAssay(obj, QFeatures::filterNA(obj[[2]],  pNA = 0), name='filtered')
#' se <- t_test_sam(obj[[3]], colData(obj), FUN = compute.t.test)
#' ind_logFC <- grep('_logFC', colnames(metadata(se)$t_test))
#' logFC <- setNames(as.data.frame(metadata(se)$t_test[,ind_logFC]), 
#' colnames(metadata(se)$t_test)[ind_logFC])
#' ind_pval <- grep('_pval', colnames(metadata(se)$t_test))
#' pval <- setNames(as.data.frame(metadata(se)$t_test[,ind_pval]), 
#' colnames(metadata(se)$t_test)[ind_pval])
#' diffAnaComputeFDR(logFC[,1], pval[,1])
#' 
#' @importFrom cp4p adjust.p
#' 
#' @export
#' 
diffAnaComputeFDR <- function(logFC, 
                              pval, 
                              threshold_PVal=0, 
                              threshold_LogFC = 0, 
                              pi0Method=1){
  
  if (missing(logFC)){
    stop("'logFC' is required.")
    return(NULL)
  }
  if (missing(pval)){
    stop("'pval' is required.")
    return(NULL)
  }
  
  upItems <- which(abs(logFC) >= threshold_LogFC)
  selectedItems <- pval[upItems]
  padj <- cp4p::adjust.p(selectedItems,  pi0Method)
  items <- which(-log10(padj$adjp[,1]) >= threshold_PVal)
  BH.fdr <- max(padj$adjp[items,2])
  
  return(BH.fdr)
}




#' @title Returns list that contains a list of the statistical tests performed with DaparToolshed and recorded
#' in an object of class \code{QFeatures}. 
#' 
#' @description This method returns a list of the statistical tests performed with xxx and formatted as a 
#' list of two DataFrame: one for the logFC data and one for the p-values data.
#' 
#' @param obj An object of class \code{SummarizedExperiment}.
#' 
#' @return A list of two slots: logFC and P_Value. Each slot contains a DataFrame 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' sTab <- colData(object)
#' gttest.se <- t_test_sam(object[[3]], sTab, FUN = compute.t.test)
#' object <- addAssay(object, gttest.se, name='t-test')
#' comp <- '25fmol_vs_10fmol'
#' object <- diff_analysis_sam(object, 't-test', name='diffAna', comp)
#' 
#' allComp <- Get_AllComparisons(object[['diffAna']])
#' 
#' 
#' @export
#' 
Get_AllComparisons <- function(obj){
  
  if(is.null(metadata(obj)$t_test)){
    stop("The  object does not contain any t-test result.")
  }
  
  logFC_KEY <- "_logFC"
  pvalue_KEY <-"_pval"
  
  ####### SAVE ALL THEPAIRWISE COMPARISON RESULTS
  res_AllPairwiseComparisons <- NULL
  ind_logFC <- grep(logFC_KEY,names(metadata(obj)$t_test))
  ind_pval <- grep(pvalue_KEY,names(metadata(obj)$t_test))
  
  #If there are already pVal values, then do no compute them 
  if (length(ind_logFC) > 0){
    res_AllPairwiseComparisons <- list(logFC = setNames(DataFrame(metadata(obj)$t_test[,ind_logFC]), names(metadata(obj)$t_test)[ind_logFC]),
                                       P_Value = setNames(DataFrame(metadata(obj)$t_test[,ind_pval]), names(metadata(obj)$t_test)[ind_pval])
    )
    
  }
  
  return(res_AllPairwiseComparisons)
}




#' @title Returns a \code{QFeatures} object with the results of
#' the differential analysis performed with 'limma package'. 
#' 
#' @description This method returns a class \code{QFeatures} object with the results
#' of differential analysis
#' 
#' @param obj An object of class \code{SUmmarizeExperiment}.
#' 
#' @param comp The name of the comparison (given by the t-tests) on whihc to apply thresholds 
#' 
#' @param th_pval xxx
#' 
#' @param th_logFC xxx
#' 
#' @return An object of class \code{SummarizedExperiment}
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' obj <- addAssay(obj, QFeatures::filterNA(obj[[2]],  pNA = 0), name='filtered')
#' obj <- addAssay(obj, t_test_sam(obj[[3]], colData(obj), FUN = 'compute.t.test'), name='t-test')
#' comp <- '25fmol_vs_10fmol'
#' da.se <- diffAnalysis(obj[['t-test']], comp, th_pval=0, th_logFC=0)
#' 
#' @export
#' 
diffAnalysis <- function(obj, 
                         comp, 
                         th_pval=0, 
                         th_logFC=0){
  
  temp <- obj
  df <- DataFrame(Significant=rep(0,nrow(obj)))
  rownames(df) <- names(obj)
  
  ##setSignificant info
  x <- metadata(obj)$t_test[,paste0(comp, '_logFC')]
  y <- -log10(metadata(obj)$t_test[,paste0(comp, '_pval')])
  
  ipval <- which(y >= th_pval)
  ilogfc <- which(abs(x) >= th_logFC)
  df$Significant[intersect(ipval, ilogfc)] <- 1
  
  
  metadata(temp)$Significant <- df
  metadata(temp)$Params$comp <- comp
  metadata(temp)$Params$th_pval <- th_pval
  metadata(temp)$Params$th_logFC <- th_logFC
  
  return(temp)
}



