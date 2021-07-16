#' @title List of hypothesis tests methods
#' 
#' @name HypothesisTestMethods
#' 
#' @export
#'
HypothesisTestMethods <- function()
  c('compute.t.test',
    'compute.group.t.test',
    'limma.complete.test',
    'wrapperClassic1wayAnova')






#' @title xxxxxx
#' 
#' @param obj xxxxx
#' 
#' @param sampleTab xxxx
#' 
#' @param logFC A vector (or list of vectors) xxxx
#' 
#' @param contrast Indicates if the test consists of the comparison of each 
#' biological condition versus 
#' each of the other ones (contrast=1; 
#' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#' 
#' @param type xxxxx
#' 
#' @return A DataFrame which contains the 
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' obj <- addAssay(obj, QFeatures::filterNA(obj[[2]],  pNA = 0), name='filtered')
#' obj <- SetAdjMat(obj, 3)
#' sTab <- colData(obj)
#' gttest <- compute.group.t.test(obj[[3]], sTab)
#' 
#' @export
#' 
#' @importFrom utils combn
#' 
compute.group.t.test <- function(obj, 
                                 sampleTab, 
                                 logFC = NULL, 
                                 contrast = "OnevsOne", 
                                 type = "Student"){
  
  if(missing(obj))
    stop("'obj' is required.")
  if(missing(sampleTab))
    stop("'sampleTab' is required.")
  if (!(contrast %in% c('OnevsOne', 'OnevsAll')))
    stop("'contrast' must be one of the following: 'OnevsOne' or 'OnevsAll")
  
  qData <- assay(obj)
  X <- GetAdjMat(obj)$all
  X.spec <- GetAdjMat(obj)$onlySpec
  
  
  .type <- type =='Student'
  sampleTab <- as.data.frame(sampleTab)
  
  res.tmp <- list()
  logFC <- list()
  P_Value <- list()
  
  sampleTab.old <- sampleTab
  Conditions.f <- factor(sampleTab$Condition, levels = unique(sampleTab$Condition))
  sampleTab <- sampleTab[unlist(lapply(split(sampleTab, Conditions.f), function(x) {x['Sample.name']})),]
  qData <- qData[,unlist(lapply(split(sampleTab.old, Conditions.f), function(x) {x['Sample.name']}))]
  Conditions <- sampleTab$Condition
  
  y <- qData
  y <- y[rownames(X), ]
  ## Keep track of proteins that are lost in aggregation step
  unsup.prot <- !(colnames(X) %in% colnames(X.spec))
  
  
  if(contrast=="OnevsOne"){
    comb <- utils::combn(levels(Conditions.f), 2)
    for(i in 1:ncol(comb)){
      c1Indice <- which(Conditions==comb[1,i])
      c2Indice <- which(Conditions==comb[2,i])
      res.tmp <- groupttest(X.spec, qData[,c1Indice], qData[,c2Indice] )
      
      #compute logFC from the result of t.test function
      p.tmp <- unlist(lapply(res.tmp,function(x) x["p.value"]))
      m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(unlist(x)["estimate.mean of x"])))
      m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(unlist(x)["estimate.mean of y"])))
      m1.name <- comb[1,i]
      m2.name <- comb[2,i]
      logFC.tmp <- m1.tmp - m2.tmp
      if (grepl(comb[1,i], m2.name)){logFC.tmp <- -logFC.tmp}
      
      peptide.spec.based.pv <- rep(1, ncol(X))
      peptide.spec.based.pv[!unsup.prot] <- p.tmp
      peptide.spec.based.FC <- rep(1, ncol(X))
      peptide.spec.based.FC[!unsup.prot] <- logFC.tmp
      
      
      txt <- paste(comb[1,i],"_vs_",comb[2,i], sep="")
      
      logFC[[paste(txt, "logFC", sep="_")]] <- peptide.spec.based.FC
      P_Value[[paste(txt, "pval", sep="_")]] <- peptide.spec.based.pv
    }
  } ##end Contrast==1
  
  if(contrast=="OnevsAll"){
    nbComp <- length(levels(Conditions.f))
    
    for(i in 1:nbComp){
      
      c1Indice <- which(Conditions==levels(Conditions.f)[i])
      # Cond.t.all <- c(1:length(Conditions))
      # Cond.t.all[c1Indice] <- levels(Conditions.f)[i]
      # Cond.t.all[-c1Indice] <- "all"
      #c1Indice <- which(Conditions==levels(Conditions.f)[i])
      res.tmp <- groupttest(X.spec, qData[,c1Indice], qData[,-c1Indice] )
      
      
      p.tmp <- unlist(lapply(res.tmp,function(x) x["p.value"]))
      m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(unlist(x)["estimate.mean of x"])))
      m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(unlist(x)["estimate.mean of y"])))
      m1.name <- levels(Conditions.f)[i]
      m2.name <-  paste("all-(",levels(Conditions.f)[i],")", sep="")
      logFC.tmp <- m1.tmp - m2.tmp
      if (grepl(levels(Conditions.f)[i], m2.name)){logFC.tmp <- -logFC.tmp}
      
      peptide.spec.based.pv <- rep(1, ncol(X))
      peptide.spec.based.pv[!unsup.prot] <- p.tmp
      peptide.spec.based.FC <- rep(1, ncol(X))
      peptide.spec.based.FC[!unsup.prot] <- logFC.tmp
      
      
      
      txt <- paste(levels(Conditions.f)[i], "_vs_(all-", levels(Conditions.f)[i], ")", sep="")
      logFC[[paste(txt, "logFC", sep="_")]] <- peptide.spec.based.FC
      P_Value[[paste(txt, "pval", sep="_")]] <- peptide.spec.based.pv
    }
  } # End Contrast=2
  
  
  res.l <- DataFrame(logFC, P_Value)
  colnames(res.l) <- c(names(logFC), names(P_Value))
  
  return(res.l) 
  
}








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
#' @details xxx
#'
#'
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept[1:1000,]
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' 
#' object <- t_test_sam(object, i=3, name = "ttest", FUN = 'compute.t.test', contrast = 'OnevsAll')
#' object <- t_test_sam(object, i=3, name = "ttest_group", FUN = 'compute.group.t.test', contrast = 'OnevsAll')
#' object <- t_test_sam(object, i=3, name = "ttest_Limma", FUN = 'limma.complete.test', comp.type= 'OnevsAll')
#' object <- t_test_sam(object, i=3, name = "ttest_Anova", FUN = 'wrapperClassic1wayAnova', with_post_hoc=TRUE, post_hoc_test='Dunnett')
#' 
"t_test_sam"

#' @param  object An object of class `SummarizedExperiment`.
#' 
#' @param sampleTab xxxxxxx
#'
#' @param FUN xxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @rdname t_test_sam
#'
setMethod("t_test_sam", "SummarizedExperiment",
          function(object,
                   sampleTab,
                   FUN,
                   ...) {
            argg <- c(as.list(environment()), list(...))
            df <- do.call(FUN, list(object, sampleTab, ...))
            
            metadata(object)$t_test <- df
            metadata(object)$Params <- argg[-match(c('object', 'sampleTab'), names(argg))]
            object
          })


#' @param  object An object of class `QFeatures`.
#' 
#' @param i A numeric vector or a character vector giving the index or the 
#'     name, respectively, of the assay(s) to be processed.
#'
#' @param name A `character(1)` naming the new assay name. Defaults
#'     are `ttestAssay`.
#' 
#' @param FUN xxx
#' 
#' @param ... Additional parameters passed to inner functions.
#' 
#' @export
#' 
#' @rdname t_test_sam
#'
setMethod("t_test_sam", "QFeatures",
          function(object, i, name = "ttestAssay", FUN,  ...) {
            if (missing(i))
              stop("Provide index or name of assay to be processed")
            if (length(i) != 1)
              stop("Only one assay to be processed at a time")
            if (is.numeric(i)) i <- names(object)[[i]]
            if (!(FUN %in% HypothesisTestMethods()))
              stop(paste0("'FUN' must be one of the following:", HypothesisTestMethods()))
            if (is.null(GetAdjMat(object[[i]])) && GetTypeDataset(object[[i]]) == 'peptide')
              object <- SetAdjMat(object, i)
            
            object <- addAssay(object,
                               t_test_sam(object[[i]], sampleTab = colData(object), FUN, ...),
                               name)
            addAssayLinkOneToOne(object, from = i, to = name)
          })






#' @title xxxxxx
#' 
#' @param X xxxxx.
#' 
#' @param qData1 xxxx.
#' 
#' @param qData2 xxxx.
#'
#' @return xxxxx
#'
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000,]
#' obj <- addAssay(obj, QFeatures::filterNA(obj[[2]],  pNA = 0), name='filtered')
#' obj <- SetAdjMat(obj, 3)
#' qData <- assay(obj[['filtered']])
#' X <- GetAdjMat(obj[[3]])$onlySpec
#' gttest <- groupttest(X, qData[,1:3], qData[,4:6])
#' 
#' @export
#' 
groupttest <- function(X, qData1=NULL, qData2 = NULL){
  res <- list()
  if (is.null(qData1) || is.null(qData2)){
    stop("At least, one condition is empty.")
  }

  for(i in 1:dim(X)[2]){
    index <- names(which(X[,i]==1))
    if(length(index)> 0){
      res[[i]] <- t.test(x=qData1[index,], y=qData2[index,], var.equal=TRUE)
    } else {
      res[[i]] <- NA
    }
  }
  #browser()
  names(res) <- colnames(X)
  return(res)
}












## Build design matrix X
# X <- as.matrix(BuildAdjacencyMatrix(obj, 'Protein_group_IDs', unique = FALSE))
# X.spec <- as.matrix(BuildAdjacencyMatrix(obj, 'Protein_group_IDs', unique = TRUE))
# 
# y <- exprs(obj)
# y <- y[rownames(X), ]
# n1 <- n2 <- 3
# n <- n1+n2 # number of samples in c1 and c2
# 
# ## Keep track of proteins that are lost in aggregation step
# unsup.prot <- !(colnames(X) %in% colnames(X.spec))
# q <- nrow(X) # Number of peptides
# p <- ncol(X) # Number of proteins
# 
# 
# # Two ways to compute it (function groupttest below)
# peptide.spec.based.tmp <- apply(X.spec, 2, FUN=function(mask) t.test(x=as.vector(y[mask == 1, 1:n1]), y=as.vector(y[mask == 1, -(1:n1)]), var.equal=TRUE)$p.value)
# peptide.spec.based.tmp <- groupttest(X.spec,y)
# 
# # then:
# peptide.spec.based.pv <- rep(1, ncol(X))
# peptide.spec.based.pv[!unsup.prot] <- peptide.spec.based.tmp
# 





# groupttest <- function(MatAdj, expr){
#   nProt <- dim(MatAdj)[2]
#   res <- rep(0,nProt)
#   
#   for(i in 1:nProt){
#     #print(i)
#     index <- names(which(MatAdj[,i]==1))
#     if(length(index)== 0){
#       res[i] <- 1
#     } else{
#       peptidestotest <- expr[index,]
#       if(length(index)== 1){
#         res[i] <- t.test(x=peptidestotest[1:3], y=peptidestotest[-(1:3)], var.equal=TRUE)$p.value
#       } else{
#         res[i] <- t.test(x=peptidestotest[,1:3], y=peptidestotest[,-(1:3)], var.equal=TRUE)$p.value
#       }
#     }
#   }
#   return(res)
# }  

