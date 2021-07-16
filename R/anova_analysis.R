
#' @title Wrapper for One-way Anova statistical test
#' 
#' @author Hélène Borges, Samuel Wieczorek
#' 
#' @param obj An object of class \code{SummarizedExperiment}.
#' 
#' @param sTab xxx
#' 
#' @param with_post_hoc a boolean saying if function must perform a Post-Hoc test or not.
#' 
#' @param post_hoc_test character string, possible values are NULL (for no
#' test; default value) or TukeyHSD" or "Dunnett". See details of
#' \code{postHocTest()} function to choose the appropriate one.
#' 
#' @details This function allows to perform a 1-way Analysis of Variance. Also
#' computes the post-hoc tests if the \code{with_post_hoc} parameter is set to
#' yes. There are two possible post-hoc tests: the Tukey Honest Significant Differences
#' (specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").
#' 
#' @return A list of two dataframes. First one called "logFC" contains
#' all pairwise comparisons logFC values (one column for one comparison) for
#' each analysed feature (Except in the case without post-hoc testing, for
#' which NAs are returned.); The second one named "P_Value" contains the
#' corresponding p-values.
#' 
#' @examples
#' library(QFeatures)
#' utils::data(Exp1_R25_prot, package='DAPARdata2')
#' obj <- Exp1_R25_prot
#' sTab <- colData(obj)
#' anovatest <- wrapperClassic1wayAnova(obj[[2]], sTab, with_post_hoc=TRUE, post_hoc_test='Dunnett')
#' 
#' @seealso [postHocTest()]
#' 
#' @export
#' 
#' @importFrom dplyr select
#' @importFrom stats aov
#' 
wrapperClassic1wayAnova <- function(obj, 
                                    sTab, 
                                    with_post_hoc = FALSE, 
                                    post_hoc_test = NULL){
 # browser()
  
  if(class(obj) != 'SummarizedExperiment')
    stop("'obj' must be of class 'SummarizedExperiment'")
  
  if (!is.null(post_hoc_test) && !(post_hoc_test %in% c('TukeyHSD', 'Dunnett')))
    stop("'post_hoc_test' must be one of the following: 'TukeyHSD', 'Dunnett'.")
  
  # if (there are NA)
  #   stop("'obj' contains missing values. Abort.")
  
  classic1wayAnova <- function(current_line, conditions){
    # vector containing the protein/peptide intensities
    intensities <- unname(unlist(current_line))
    # anova sur ces deux vecteurs
    aov_test <- aov(formula = intensities ~ conditions, data = NULL)
    return(aov_test)
    
  }
  
  
  qData <- as.data.frame(assay(obj))
  
  if (!(with_post_hoc %in% c(TRUE, FALSE)))
    stop("Wrong with_post_hoc parameter. Please choose between FALSE or TRUE")
  
  if (!isTRUE(with_post_hoc) && is.null(post_hoc_test))
    stop("You want to perform a post-hoc test but did not specify which test. Please choose between 'Dunnett' or 'TukeyHSD'.")
  
   
  if (!isTRUE(with_post_hoc)){
    anova_tests <- as.data.frame(t(apply(qData, 1, function(x) unlist(summary(classic1wayAnova(x, conditions=as.factor(sTab$Condition)))))))
    results <- dplyr::select(anova_tests, `Pr(>F)1`)
    to_return <- DataFrame("logFC" = data.frame("anova_1way_logFC" = matrix(NA, nrow = nrow(results)), row.names = rownames(results)),
                      "P_Value" = data.frame("anova_1way_pval" = results$`Pr(>F)1`, row.names = rownames(results)))
  } else {
      anova_tests <- t(apply(qData, 1, classic1wayAnova, conditions = as.factor(sTab$Condition)))
      names(anova_tests) <- rownames(qData)
      to_return <- postHocTest(aov_fits = anova_tests,
                               post_hoc_test = post_hoc_test)
  } 
  
  return(to_return)
}



#' Extract logFC and raw pvalues from multiple post-hoc models summaries
#'
#' @param post_hoc_models_summaries a list of summaries of post-hoc models.
#'
#' @return a list of 2 dataframes containing the logFC values and pvalues for
#' each comparison.
#' 
#' @author Hélène Borges
#' 
#' @examples
#' 
#' @importFrom purrr map_dfr
#' 
formatPHResults <- function(post_hoc_models_summaries){
  
  # récupérer les différences entre les moyennes
  res_coeffs <- lapply(post_hoc_models_summaries, function(x) x$test$coefficients)
  logFC <- data.frame(purrr::map_dfr(res_coeffs, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients))
  logFC <- as.data.frame(t(logFC))
  # extract raw p-values (non-adjusted)
  res_pvals <- lapply(post_hoc_models_summaries, function(x) x$test$pvalues)
  pvals <- data.frame(purrr::map_dfr(res_pvals, cbind),
                      row.names = names(post_hoc_models_summaries[[1]]$test$coefficients))
  pvals <- as.data.frame(t(pvals))
  
  res <- DataFrame(logFC, pvals)
  names(logFC) <- stringr::str_c(stringr::str_replace(colnames(logFC), " - ", "_vs_"), "_logFC")
  names(pvals) <- stringr::str_c(stringr::str_replace(colnames(pvals), " - ", "_vs_"), "_pval")
  colnames(res) <- c(names(logFC), names(pvals))
  
  return(res)
}





#' Post-hoc tests for classic 1-way ANOVA
#'
#' @description This function allows to compute a post-hoc test after a 1-way
#' ANOVA analysis. It expects as input an object obtained with the function
#' \code{classic1wayAnova}. The second parameter allows to choose between 2
#' different post-hoc tests: the Tukey Honest Significant Differences
#' (specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").
#' 
#'
#' @param aov_fits a list containing aov fitted model objects
#' 
#' @param post_hoc_test a character string indicating which post-hoc test to
#' use. Possible values are "TukeyHSD" or "Dunnett". See details for what to
#' choose according to your experimental design.
#' 
#' @details
#' This is a function allowing to realise post-hoc tests for a set of
#' proteins/peptides for which a classic 1-way anova has been performed with
#' the function \code{classic1wayAnova}. Two types of tests are currently
#' available: The Tukey HSD's test and the Dunnett's test. Default is Tukey's
#' test.
#' The Tukey HSD's test compares all possible pairs of means, and is based on a
#' studentized range distribution. Here is used the \code{TukeyHSD()} function,
#' which can be applied to balanced designs (same number of samples in each
#' group), but also to midly unbalanced designs.
#' The Dunnett's test compares a single control group to all other groups.
#' Make sure the factor levels are properly ordered.
#' 
#' @return a list of 2 dataframes: first one called "LogFC" contains
#' all pairwise comparisons logFC values (one column for one comparison) for
#' each analysed feature; The second one named "P_Value" contains the
#' corresponding pvalues.
#' 
#' @author Hélène Borges
#' 
#' @examples 
#' 
#' 
#' 
#' @export
#' 
#' @importFrom multcomp glht adjusted mcp
#' 
postHocTest <- function(aov_fits, post_hoc_test = "TukeyHSD"){
  
  # Check post_hoc_test parameter
  if (!(post_hoc_test %in% c('TukeyHSD', 'Dunnett')))
  stop("Wrong post_hoc_test parameter. Please choose between TukeyHSD or
             Dunnett.")
             
  if (post_hoc_test == "TukeyHSD")
    post_hoc_test <- 'Tukey'
  
  # use of adjusted("none") to obtain raw p-values (and not adjusted ones)
    models_summaries <- lapply(aov_fits,
                               function(x) summary(multcomp::glht(x, linfct = multcomp::mcp(conditions = post_hoc_test)),
                                                       test = multcomp::adjusted("none")))
   
    res <- formatPHResults(models_summaries)
  
  return(res)
}