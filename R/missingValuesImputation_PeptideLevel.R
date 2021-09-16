


#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param x xxxxx
#' 
#' @param sampleTab xxxx
#' 
#' @return The \code{assay(obj)} matrix with imputed values instead of missing values.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' sampleTab <- colData(obj)
#' dat <- impute_mle_dapar(assay(obj[[2]]), colData(obj))
#' 
#' @export
#' 
#' @importFrom imp4p impute.mle
#' 
impute_mle_dapar <- function(x, sampleTab){
  
  if (missing(sampleTab))
    stop("'sampleTab' is missing.")
  
  # sort conditions to be compliant with impute.mle (from imp4p version 0.9) which only manage ordered samples 
  old.sample.name <- sampleTab$Sample.name
  new.order <- order(sampleTab$Condition)
  new.sampleTab <- sampleTab[new.order,]
  conds <- factor(new.sampleTab$Condition, levels=unique(new.sampleTab$Condition))
  new.x <- x[,new.order]
  res <- imp4p::impute.mle(as.matrix(new.x), conditions=conds)
  
  #restore old order
  res <- res[,old.sample.name]
  return (res)
}





#' @title Missing values imputation using the LSimpute algorithm.
#' 
#' @param x xxxxx
#' 
#' @param sampleTab xxxx
#' 
#' @param lapala xxxxxxxxxxx
#' 
#' @param progress.bar Same as the function \code{mi.mix} in the package \code{imp4p}
#' 
#' @param ... xxxxxxxx
#' 
#' @return The \code{assay(obj)} matrix with imputed values instead of missing values.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[seq_len(1000),]
#' x <- assay(obj[[2]])
#' sampleTab <- colData(obj)
#' dat <- impute_mi(x, sampleTab)
#' 
#' @export
#' 
#' @importFrom imp4p impute.rand prob.mcar.tab estim.mix estim.bound
#' 
impute_mi <- function (x, sampleTab, lapala = TRUE, progress.bar=TRUE, ...) 
{
  
  if (missing(sampleTab))
    stop("'sampleTab' is missing.")
  
  # sort conditions to be compliant with impute.mle (from imp4p version 0.9) which only manage ordered samples 
  old.sample.name <- sampleTab$Sample.name
  new.order <- order(sampleTab$Condition)
  new.sampleTab <- sampleTab[new.order,]
  conds <- factor(new.sampleTab$Condition, levels=unique(new.sampleTab$Condition))
  new.x <- x[,new.order]
  
  conditions <- as.factor(new.sampleTab$Condition)
  repbio <- as.factor(new.sampleTab$Bio.Rep)
  reptech <-as.factor(new.sampleTab$Tech.Rep)
  
  
  tab <- new.x
  
  if (progress.bar == TRUE) {
    cat(paste("\n 1/ Initial imputation under the MCAR assumption with impute.rand ... \n  "))
  }
  dat.slsa = imp4p::impute.rand(tab = tab, conditions = conditions)
  
  if (progress.bar == TRUE) {
    cat(paste("\n 2/ Estimation of the mixture model in each sample... \n  "))
  }
  res = imp4p::estim.mix(tab = tab, 
                         tab.imp = dat.slsa, 
                         conditions = conditions, 
                         x.step.mod = 300, 
                         x.step.pi = 300, 
                         nb.rei = 100)
  
  
  if (progress.bar == TRUE) {
    cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "))
  }
  born = imp4p::estim.bound(tab = tab, 
                            conditions = conditions, 
                            q = 0.95)
  
  proba = imp4p::prob.mcar.tab(born$tab.upper, res)
  
  
  if (progress.bar == TRUE) {
    cat(paste("\n 4/ Multiple imputation strategy with mi.mix ... \n  "))
  }
  data.mi = imp4p::mi.mix(tab = tab, 
                          tab.imp = dat.slsa, 
                          prob.MCAR = proba, 
                          conditions = conditions, 
                          repbio = repbio, 
                          reptech = reptech, 
                          ...)
  
  colnames(data.mi) <- new.sampleTab$Sample.name
  rownames(data.mi) <- rownames(new.x)  
  
  if (lapala == TRUE){
    if (progress.bar == TRUE) {
      cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa ... \n  "))
    }
    data.final <- impute_pa2(data.mi, 
                             new.sampleTab, 
                             ...)
  } else {
    data.final <- data.mi
  }
  
  # restore previous order
  data.final <- data.final[,old.sample.name]
  return (data.final)
  
}


#' @title Generator of simulated values
#' 
#' @param n An integer which is the number of simulation (same as in rbeta)
#' 
#' @param min An integer that corresponds to the lower bound of the interval
#' 
#' @param max An integer that corresponds to the upper bound of the interval
#' 
#' @param param1 An integer that is the first parameter of rbeta function.
#' 
#' @param param2 An integer that is second parameter of rbeta function.
#' 
#' @return A vector of n simulated values
#' 
#' @author Thomas Burger
#' 
#' @examples
#' translatedRandomBeta(1000, 5, 10, 1, 1)
#' 
#' @export
#' 
#' @importFrom stats rbeta
#' 
translatedRandomBeta <- function(n, min, max, param1=3, param2=1){
  if(missing(n))
    stop("'n' is missing.")
  if(missing(min))
    stop("'min' is missing.")
  if(missing(max))
    stop("'max' is missing.")
  
  
  scale <- max - min
  simu <- stats::rbeta(n, param1, param2)
  res <- (simu*scale) + min
  return(res)
}









################################################
#' This method is a variation to the function \code{impute.pa} from the package 
#' \code{imp4p}.
#' 
#' @title Missing values imputation from a \code{MSnSet} object
#' 
#' @param x xxxxx
#' 
#' @param sampleTab xxxxx
#' 
#' @param q.min A quantile value of the observed values allowing defining the 
#' maximal value which can be generated. This maximal value is defined by the
#' quantile q.min of the observed values distribution minus eps. 
#' Default is 0 (the maximal value is the minimum of observed values minus eps).
#' 
#' @param q.norm A quantile value of a normal distribution allowing defining 
#' the minimal value which can be generated. Default is 3 (the minimal value 
#' is the maximal value minus qn*median(sd(observed values)) where sd is the 
#' standard deviation of a row in a condition).
#' 
#' @param eps A value allowing defining the maximal value which can be 
#' generated. This maximal value is defined by the quantile q.min of the 
#' observed values distribution minus eps. Default is 0.
#' 
#' @param distribution The type of distribution used. Values are unif or beta.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[seq_len(1000)]
#' x <- assay(obj[[2]])
#' sampleTab <- colData(obj)
#' dat <- impute_pa2(x, sampleTab,distribution="beta")
#' 
#' @export
#' 
#' @importFrom stats median quantile runif sd
#' 
#' @return NA
#' 
impute_pa2 <- function (x, sampleTab, q.min = 0, q.norm = 3, eps = 0, distribution = "unif"){
  
  if (missing(sampleTab))
    stop("'sampleTab' is missing.")
  
  if (missing(x))
    stop("'x' is missing.")
  
  # sort conditions to be compliant with impute.mle (from imp4p version 0.9) which only manage ordered samples     old.sample.name <- sampleTab$Sample.name
  old.sample.name <- sampleTab$Sample.name
  new.order <- order(sampleTab$Condition)
  new.sampleTab <- sampleTab[new.order,]
  conds <- factor(new.sampleTab$Condition, levels=unique(new.sampleTab$Condition))
  tab_imp <- x[,new.order]
  
  conditions <-  as.factor(new.sampleTab$Condition)
  
  qu = apply(tab_imp, 2, quantile, na.rm = TRUE, q.min)
  nb_cond = length(levels(conditions))
  nb_rep = rep(0, nb_cond)
  k = 1
  j = 1
  for (i in seq_len(nb_cond)) {
    nb_rep[i] = sum((conditions == levels(conditions)[i]))
    sde = apply(tab_imp[, (k:(k + nb_rep[i] - 1))], 1, sd, 
                na.rm = TRUE)
    while (j < (k + nb_rep[i])) {
      if(distribution == "unif")
      {
        tab_imp[which(is.na(tab_imp[, j])), j] = stats::runif(n = sum(is.na(tab_imp[,j])), 
                                                              min = qu[j] - eps - q.norm * median(sde, na.rm = TRUE), 
                                                              max = qu[j] - eps)
      } else if (distribution == "beta"){
        tab_imp[which(is.na(tab_imp[, j])), j] = translatedRandomBeta(n = sum(is.na(tab_imp[,j])), 
                                                                      min = qu[j] - eps - q.norm * median(sde, na.rm = TRUE), 
                                                                      max = qu[j] - eps,
                                                                      param1 = 3,
                                                                      param2 = 1)
      }
      j = j + 1
    }
    k = k + nb_rep[i]
  }
  
  tab_imp <- tab_imp[,old.sample.name]
  return(tab_imp)
}










