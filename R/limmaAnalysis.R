#' @title Check the validity of the experimental design
#' 
#' @param sTab xxxx
#' 
#' @return A boolean
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' CheckDesign(colData(Exp1_R25_pept))
#' 
#' @export
#' 
CheckDesign <- function(sTab){
  
  res <- list(valid = FALSE,
              warn = NULL)
  
  names <- colnames(sTab)
  level.design <- sum(names %in% c('Bio.Rep', 'Tech.Rep', 'Analyt.Rep'))
  
  
  res <- check.conditions(sTab$Condition)
  if (!res$valid){
    return(res)
  }
  # Check if all the column are fullfilled
  
  if (level.design == 1){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE, warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    }
  }
  else if (level.design == 2){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE, warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    }else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)){
      res <- list(valid=FALSE, warn="The Tech.Rep colmumn are not full filled.")
      return(res)
    }
  }
  else if (level.design == 3){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE, warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)){
      res <- list(valid=FALSE, warn="The Tech.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Analyt.Rep) || (NA %in% sTab$Analyt.Rep)){
      res <- list(valid=FALSE, warn="The Analyt.Rep colmumn are not full filled.")
      return(res)
    }
  }
  
  # Check if the hierarchy of the design is correct
  if (level.design == 1){res <- test.design(sTab[,c("Condition", "Bio.Rep")])}
  else if (level.design == 2){res <- test.design(sTab[,c("Condition", "Bio.Rep","Tech.Rep")])}
  else if (level.design == 3){
    res <- test.design(sTab[,c("Condition", "Bio.Rep","Tech.Rep")])
    if (res$valid)
    {
      res <- test.design(sTab[,c("Bio.Rep","Tech.Rep", "Analyt.Rep")])
      
    }
  }
  
  return(res)
  
  
  }




#' @title Check if xxxxxx
#' 
#' @param tab A DataFrame which correspond to the samples definition. This DataFrame contains three columns: 
#' "Bio.Rep" , "Tech.Rep" and "Analyt.Rep"
#' 
#' @return A list of two items:
#' * xxxxxx
#' * xxxxxxx
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' test.design(colData(Exp1_R25_pept)[,1:3])
#' 
#' @export
#' 
test.design <- function(tab){
  valid <- TRUE
  txt <- NULL
  level <- NULL
  
  level.a  <- factor(tab[,1], ordered = TRUE)
  level.b <- factor(tab[,2], ordered = TRUE)
  name.level.a <- colnames(tab)[1]
  name.level.b <-colnames(tab)[2]
  
  level.c <- NULL
  level.c <- if(ncol(tab)==3) factor(tab[,3], ordered = TRUE) else NULL
  name.level.c <- if(ncol(tab)==3) colnames(tab)[3] else NULL
  
  
  # verification intersection sur B
  ##verification de la non redondance'intersection vide entre les groupes
  uniqueA <- unique(level.a)
  ll <- lapply(uniqueA, function(x){as.character(level.b)[which(level.a==x)]})
  n <- NULL
  for (i in 1:(length(uniqueA)-1)){
    for (j in (i+1):length(uniqueA)){
      n <- c(n,intersect(ll[[i]], ll[[j]]))
    }
  }
  if (length(n) > 0){
    valid <- FALSE
    txt <- c(txt,paste0("The value ", n, " in column '", colnames(tab)[2], "' is not correctly set.\n"))
  }
  
  
  #verification si niveau hierarchique inf
  if (length(levels(level.a)) == length(levels(level.b))){
    ## c'est un design de niveau n-1 en fait
    valid <- FALSE
    txt <- c(txt,paste0("The column ",name.level.b, " is not informative. Thus, the design is not of level (n-1).\n"))
  } 
  else if (!is.null(level.c)){
    if (length(levels(level.b)) == length(levels(level.c))){
      ## c'est un design de niveau n-1 en fait
      valid <- FALSE
      txt <- c(txt,paste0("The column ",name.level.c, " is not informative. Thus, the design is of level (n-1).\n"))
    }
  } 
  
  #verification si niveau non informatif
  return(list(valid = valid, 
              warn = txt))
}



#' @title Check if the design of conditions is valid
#' 
#' @param conds A vector of conditions of the samples.
#' 
#' @return A list
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' check.conditions(colData(Exp1_R25_pept)$Condition)
#' 
#' @export
#' 
check.conditions <- function(conds){
  res <- list(valid=TRUE,warn=NULL)
  
  if (("" %in% conds) || (NA %in% conds)){
    res <- list(valid=FALSE,warn="The conditions are not full filled.")
    return(res)
  }
  
  # Check if there is at least two conditions
  if (length(unique(conds)) < 2){
    res <- list(valid=FALSE, warn="The design must contain at least two conditions.")
    return(res)
  }
  
  
  # check if each condition has at least two values
  nValPerCond <- unlist(lapply(unique(conds), function(x){length(conds[which(conds==x)])}))
  if (all(nValPerCond < 2)){
    res <- list(valid=FALSE,warn="The design must contain at least two values per condition.")
    return(res)
  }
  
  return(res)
}


#' @title Builds the design matrix
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' make.design(colData(Exp1_R25_pept))
#' 
#' @export
#' 
make.design <- function(sTab){
  
  if (!CheckDesign(sTab)$valid){
    warning("The design matrix is not correct.")
    warning(CheckDesign(sTab)$warn)
    return(NULL)
  }
  
  n <- ncol(sTab)
  if (n <= 2){
    stop("Error in design matrix dimensions which must have at least 3 columns.")
  } 
  
  res <- do.call(paste0("make.design.", (n-2)),list(sTab))
  
  return(res)
}





#' @title Builds the design matrix for designs of level 1
#' 
#' @param sTab The data.frame which correspond to the definition of replicates
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' make.design.1(colData(Exp1_R25_pept))
#' 
#' @export
#' 
make.design.1 <- function(sTab){
  
  Conditions <- factor(sTab$Condition, 
                       levels=unique(sTab$Condition)
                       , ordered = TRUE)
  nb_cond <- length(unique(Conditions))
  nb_samples <- nrow(sTab)
  
  #CGet the number of replicates per condition
  nb_Rep <- rep(0, nb_cond)
  for (i in 1:nb_cond){
    nb_Rep[i] = sum((Conditions==unique(Conditions)[i]))
  }
  
  design <- matrix(0, nb_samples, nb_cond)
  n0=1
  coln=NULL
  for (j in 1:nb_cond){
    coln=c(coln, paste("Condition", j, collapse=NULL,sep=""))
    design[(n0:(n0+nb_Rep[j]-1)),j] <- rep(1,length((n0:(n0+nb_Rep[j]-1))))
    n0 <- n0 + nb_Rep[j]
  }
  colnames(design) <- coln
  
  return(design)
}





#' @title Builds the design matrix for designs of level 2
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' \donttest{
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' make.design.2(colData(Exp1_R25_pept))
#' }
#' 
#' @export
#' 
#' @importFrom stats rnorm
#' 
make.design.2 <- function(sTab){
  Condition <- factor(sTab$Condition,
                      levels=unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep, 
                   levels=unique(sTab$Bio.Rep))
  
  #Renome the levels of factor
  levels(Condition) <- c(1:length(levels(Condition)))
  levels(RepBio) <- c(1:length(levels(RepBio)))
  
  #Initial design matrix
  df <- rep(0,nrow(sTab))
  names(df) <- rownames(sTab)
  design <- stats::model.matrix(df~0+Condition:RepBio)
  
  #Remove empty columns in the design matrix
  design <- design[,(apply(design,2,sum)>0)]
  #Remove identical columns in the design matrix
  coldel <- -1
  for (i in 1:(length(design[1,])-1)){
    d2 <- as.matrix(design[,(i+1):length(design[1,])]);
    for (j in 1:length(d2[1,])){
      d2[,j]=d2[,j]-design[,i];
    }
    e <- as.matrix(stats::rnorm(length(design[,1]),10,1));
    sd2 <- t(e)%*%d2
    liste <- which(sd2==0)
    coldel <- c(coldel,liste+i)
  }
  design <- design[,(1:length(design[1,]))!=coldel]
  colnames(design) <- make.names(colnames(design))
  return(design)
}


#' @title Builds the design matrix for designs of level 3
#' 
#' @param sTab The data.frame which correspond to the definitino of the replicates
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' sTab <-cbind(colData(Exp1_R25_pept), Tech.Rep=1:6)
#' make.design.3(sTab)
#' 
#' @export
#' 
#' @importFrom stats model.matrix rnorm
#' 
make.design.3 <- function(sTab){
  
  Condition <- factor(sTab$Condition,  levels=unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep,  levels=unique(sTab$Bio.Rep))
  RepTech <- factor(sTab$Tech.Rep,  levels=unique(sTab$Tech.Rep))
  
  
  #Rename the levels of factor
  levels(Condition) <- c(1:length(levels(Condition)))
  levels(RepBio) <- c(1:length(levels(RepBio)))
  levels(RepTech) <- c(1:length(levels(RepTech)))
  
  
  #Initial design matrix
  df <- rep(0,nrow(sTab))
  names(df) <- rownames(sTab)
  design <- stats::model.matrix(df~0+Condition:RepBio:RepTech)
  
  #Remove empty columns in the design matrix
  design=design[,(apply(design,2,sum)>0)]
  
  #Remove identical columns in the design matrix
  coldel <- -1
  for (i in 1:(length(design[1,])-1)){
    d2=as.matrix(design[,(i+1):length(design[1,])]);
    for (j in 1:length(d2[1,])){
      d2[,j]=d2[,j]-design[,i];
    }
    e=as.matrix(stats::rnorm(length(design[,1]),10,1));
    sd2=t(e)%*%d2
    liste=which(sd2==0)
    coldel=c(coldel,liste+i)
  }
  design <- design[,(1:length(design[1,]))!=coldel]
  colnames(design) <- make.names(colnames(design))
  return(design)
}






#' @title Builds the contrast matrix
#' 
#' @param design The data.frame which correspond to the xxxxx
#' 
#' @param condition xxxxx
#' 
#' @param contrast An integer that indicates if the test consists of the comparison of each 
#' biological condition versus each of the other ones (Contrast=1; 
#' for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#' or each condition versus all others (Contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#'  H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#'  
#' @return A constrat matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' design <- make.design(colData(Exp1_R25_pept))
#' conds <- colData(Exp1_R25_pept)$Condition
#' make.contrast(design, conds)
#' 
#' @export
#' 
make.contrast <- function(design, condition, contrast = 1){
  
  #######################################################
  aggreg.column.design <- function(design,Condition){
    nb.cond <- length(levels(Condition))
    name.col <- colnames(design)
    name.cond <- NULL
    nb.col <- NULL
    for (i in 1:nb.cond){
      col.select <- NULL
      col.name.begin <- paste("Condition",i, sep = "")
      nc <- nchar(col.name.begin)
      for (j in 1:length(design[1,])){
        if (substr(name.col[j], 1, nc) == col.name.begin){
          col.select <- c(col.select,j)
        }
      }
      name.aggreg <- NULL
      for (j in 1:length(col.select)){
        name.aggreg  <- paste(name.aggreg,name.col[col.select[j]], sep="+")
      }
      name.aggreg <- substr(name.aggreg, 2, nchar(name.aggreg))
      name.cond <- c(name.cond,name.aggreg)
      nb.col <- c(nb.col,length(col.select))
    }
    return(list(name.cond, nb.col))
  }
  
  
  
  nb.cond=length(levels(condition))
  r=aggreg.column.design(design,condition)
  label.agg=r[[1]]
  nb.agg=r[[2]]
  k=1
  
  if (contrast == 1){
    ## Contrast for One vs One
    contra=rep(0,sum(1:(nb.cond-1)))
    for (i in 1:(nb.cond-1)){
      for (j in (i+1):nb.cond){
        contra[k]=c(paste("(",label.agg[i],")/",
                          nb.agg[i],"-(",label.agg[j],")/",
                          nb.agg[j]))
        k=k+1
      }
    }
  } else if (contrast==2){
    ## Contrast for One vs All
    contra=rep(0,nb.cond)
    for (i in 1:(nb.cond)){
      contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i]))
      nb=sum(nb.agg[(1:nb.cond)[(1:nb.cond)!=i]])
      for (j in (1:nb.cond)[(1:nb.cond)!=i]){
        contra[k]=c(paste(contra[k],"-(",label.agg[j],")/",nb))
      }
      k=k+1
    }
  }
  
  return(contra)
}



#' @title Computes a hierarchical differential analysis
#' 
#' @param obj A matrix of quantitative data, without any missing values.
#' 
#' @param sTab A DataFrame of experimental design. 
#' 
#' @param comp.type A string that corresponds to the type of comparison. 
#' Values are: 'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.
#' 
#' @return A list of two dataframes : logFC and P_Value. The first one contains
#' the logFC values of all the comparisons (one column for one comparison), the second one contains
#' the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes
#' are identical and correspond to the description of the comparison. 
#' 
#' @author Hélène Borges, Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' Exp1_R25_pept <- addAssay(Exp1_R25_pept, QFeatures::filterNA(Exp1_R25_pept[[2]],  pNA = 0), name='filtered')
#' sTab <- colData(Exp1_R25_pept)
#' se <- Exp1_R25_pept[['filtered']]
#' limma <- limmaCompleteTest(se, sTab, comp.type='anova1way')
#' 
#' @importFrom limma contrasts.fit makeContrasts lmFit
#' @importFrom tidyr crossing %>%
#' @importFrom dplyr mutate filter select
#' 
#' @export
#' 
limmaCompleteTest <- function(obj,
                              sTab, 
                              comp.type = "OnevsOne"){
  
  if(class(obj) != 'SummarizedExperiment')
    stop("'obj' must be of class 'SummarizedExperiment'")
  
  if (!(comp.type %in% c('OnevsOne', 'OnevsAll', 'anova1way')))
    stop("'comp.type' must be one of the following: 'OnevsOne', 'OnevsAll', 'anova1way.")
  
  
  qData <- assay(obj)
  switch(comp.type,
         OnevsOne = contrast <- 1,
         OnevsAll = contrast <- 2)
  
  sTab <- as.data.frame(sTab)
  sTab.old <- sTab
  conds <- factor(sTab$Condition, levels=unique(sTab$Condition))
  sTab <- sTab[unlist(lapply(split(sTab, conds), function(x) {x['Sample.name']})),]
  qData <- qData[,unlist(lapply(split(sTab.old, conds), function(x) {x['Sample.name']}))]
  conds <- conds[order(conds)]
  
  
  res.l <- NULL
  design.matrix <- make.design(sTab)
  if(!is.null(design.matrix)) {
    if(comp.type == 'OnevsOne' || comp.type == "OnevsAll"){
      contra <- make.contrast(design.matrix,condition=conds, contrast)
      cmtx <- limma::makeContrasts(contrasts=contra,
                                   levels=make.names(colnames(design.matrix)))
      fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(qData, design.matrix), cmtx))
      res.l <- formatLimmaResult(fit, conds, contrast)
    }else if(comp.type == "anova1way"){
      # make the orthogonal contrasts
      contrasts <- tidyr::crossing(A = colnames(design.matrix), B = colnames(design.matrix), .name_repair = "minimal") %>%
        dplyr::filter(A!=B)
      orthogonal_contrasts <- dplyr::filter(contrasts, !duplicated(paste0(pmax(A, B), pmin(A, B)))) %>%
        dplyr::mutate(contrasts = stringr::str_glue("{A}-{B}"))
      # make the contrasts in a format adapted for limma functions
      contrasts_limma_format <- limma::makeContrasts(contrasts = orthogonal_contrasts$contrasts,
                                                     levels = colnames(design.matrix))
      ebayes_fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(qData, design.matrix), contrasts_limma_format))
      fit_table <- limma::topTable(ebayes_fit, sort.by = "none", number = nrow(qData))
      fit_pvalue <- dplyr::select(fit_table, "anova_1way_pval" = P.Value)
      res.l <- list("logFC" = data.frame("anova_1way_logFC" = matrix(NA, nrow = nrow(fit_pvalue)),
                                         row.names = rownames(fit_pvalue)),
                    "P_Value" = fit_pvalue)
    }
  }
  
  return(res.l)
}





#' @title xxxx
#' 
#' @param fit xxxx
#' 
#' @param conds xxxx
#' 
#' @param contrast xxxx
#' 
#' @return A list of two dataframes : logFC and P_Value. The first one contains
#' the logFC values of all the comparisons (one column for one comparison), the second one contains
#' the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes
#' are identical and correspond to the description of the comparison. 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' object <- Exp1_R25_pept
#' object <- addAssay(object, QFeatures::filterNA(object[[2]],  pNA = 0), name='filtered')
#' qData <- assay(object, 'filtered')
#' sTab <- colData(object)
#' limma <- limmaCompleteTest(qData, sTab)
#' }
#' 
#' @importFrom stringr str_match_all
#' 
formatLimmaResult <- function(fit, conds, contrast){
  
  #res.tmp <- topTable(fit,number=Inf, sort.by="none")
  #res <- cbind(res.tmp[,1:Compa.Nb], fit$p.value)
  #names(res) <- gsub(".", "_", names(res), fixed=TRUE)
  res <- cbind(fit$coefficients, fit$p.value)
  
  #how many comparisons have been done (and thus how many columns of pval)
  Compa.Nb <- dim(fit$p.value)[2]
  #empty colnames vector
  cn <- c()
  for (i in 1:Compa.Nb){
    
    #not the same syntax to pars if Contast=1 or Contrast=2
    if(contrast==1){
      compa <- as.data.frame(stringr::str_match_all(colnames(fit$p.value)[i],"[[:space:]]Condition([[:digit:]]+)"))
      cn[i] <- paste(unique(conds)[as.numeric(compa[1,2])], "_vs_",unique(conds)[as.numeric(compa[2,2])], sep="")
    }
    if(contrast==2){
      #hierarchic only
      #compa<-str_match_all(colnames(fit$p.value)[i], "[[:space:]]Condition([[:digit:]]+)[[:space:]]")[[1]]
      #cn[i]<-paste(levels(Conditions)[as.numeric(compa[1,2])], "vs(all-",levels(Conditions)[as.numeric(compa[1,2])], ")", sep="")
      
      #hier and non hier
      compa <- stringr::str_match_all(colnames(fit$p.value)[i], "[[:space:]]Condition([[:digit:]]+)")[[1]]
      cn[i] <- paste(unique(conds)[as.numeric(compa[1,2])], "_vs_(all-",unique(conds)[as.numeric(compa[1,2])], ")", sep="")
    }
  }
  
  # res.l <- list(
  #   logFC = as.data.frame(res[,1:Compa.Nb]),
  #   P_Value = as.data.frame(res[,-(1:Compa.Nb)] )
  # )
  # 
  # colnames(res.l$logFC) <- paste(cn, "logFC",sep="_")
  # colnames(res.l$P_Value) <- paste(cn, "pval",sep="_")
  
  res.l <- DataFrame(data.frame(res[,1:Compa.Nb], res[,-(1:Compa.Nb)]))
  colnames(res.l) <- c(paste(cn, "logFC",sep="_"), paste(cn, "pval",sep="_"))
  
  ## end colnames
  
  return(res.l)
}

