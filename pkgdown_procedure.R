setwd("Github/AdaptedForFeatures/DAPAR2_pkgdown/")
library(pkgdown)
# copy logo.png in "."

#############################################################################################
init_site()
# One time
# create docs/
# Minor, create pkgdown/favicon and 404.html

# -- Initialising site -----------------------------------------------
# Copying '../../../R/win-library/4.0/pkgdown/assets/bootstrap-toc.css' to 'bootstrap-toc.css'
# Copying '../../../R/win-library/4.0/pkgdown/assets/bootstrap-toc.js' to 'bootstrap-toc.js'
# Copying '../../../R/win-library/4.0/pkgdown/assets/docsearch.css' to 'docsearch.css'
# Copying '../../../R/win-library/4.0/pkgdown/assets/docsearch.js' to 'docsearch.js'
# Copying '../../../R/win-library/4.0/pkgdown/assets/link.svg' to 'link.svg'
# Copying '../../../R/win-library/4.0/pkgdown/assets/pkgdown.css' to 'pkgdown.css'
# Copying '../../../R/win-library/4.0/pkgdown/assets/pkgdown.js' to 'pkgdown.js'
# -- Building favicons -----------------------------------------------
# Building favicons with realfavicongenerator.net...
# Copying 'pkgdown/favicon/apple-touch-icon-120x120.png' to 'apple-touch-icon-120x120.png'
# Copying 'pkgdown/favicon/apple-touch-icon-152x152.png' to 'apple-touch-icon-152x152.png'
# Copying 'pkgdown/favicon/apple-touch-icon-180x180.png' to 'apple-touch-icon-180x180.png'
# Copying 'pkgdown/favicon/apple-touch-icon-60x60.png' to 'apple-touch-icon-60x60.png'
# Copying 'pkgdown/favicon/apple-touch-icon-76x76.png' to 'apple-touch-icon-76x76.png'
# Copying 'pkgdown/favicon/apple-touch-icon.png' to 'apple-touch-icon.png'
# Copying 'pkgdown/favicon/favicon-16x16.png' to 'favicon-16x16.png'
# Copying 'pkgdown/favicon/favicon-32x32.png' to 'favicon-32x32.png'
# Copying 'pkgdown/favicon/favicon.ico' to 'favicon.ico'
# Copying 'logo.png' to 'logo.png'
# Writing '404.html'


#############################################################################################
build_home()
# create index.html from more or less README.md (doesn't work directly with .Rmd, have to be "knit" into .md)
# change YAML for details

# -- Building home --------------------------------------------------------------------------------------------
#   Writing 'authors.html'
# -- Previewing site ------------------------------------------------------------------------------------------


#############################################################################################
# build_reference_index()
# # Writing 'reference/index.html'
build_reference()
# -- Building function reference ------------------------------------------------------------------------------
# Loading DAPAR2
# Registered S3 method overwritten by 'quantmod':
#   method            from
# as.zoo.data.frame zoo 
# 
# Loading required package: Biobase
# Loading required package: BiocGenerics
# Loading required package: parallel
# 
# Attaching package: ‘BiocGenerics’
# 
# The following objects are masked from ‘package:parallel’:
#   
#   clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply,
# parCapply, parLapply, parLapplyLB, parRapply, parSapply, parSapplyLB
# 
# The following objects are masked from ‘package:stats’:
#   
#   IQR, mad, sd, var, xtabs
# 
# The following objects are masked from ‘package:base’:
#   
#   anyDuplicated, append, as.data.frame, basename, cbind, colnames, dirname, do.call,
# duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
# mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
# Reduce, rownames, sapply, setdiff, sort, table, tapply, union, unique, unsplit, which.max,
# which.min
# 
# Welcome to Bioconductor
# 
# Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor,
# see 'citation("Biobase")', and for packages 'citation("pkgname")'.
# 
# 
# Attaching package: ‘DynDoc’
# 
# The following object is masked from ‘package:BiocGenerics’:
#   
#   path
# 
# Reading 'man/addConnexComp.Rd'
# Writing 'reference/addConnexComp.html'
# Reading 'man/addListAdjacencyMatrices.Rd'
# Writing 'reference/addListAdjacencyMatrices.html'
# Reading 'man/addOriginOfValues.Rd'
# Writing 'reference/addOriginOfValues.html'
# Reading 'man/aggIter.Rd'
# Writing 'reference/aggIter.html'
# Reading 'man/aggIterParallel.Rd'
# Writing 'reference/aggIterParallel.html'
# Reading 'man/aggMean.Rd'
# Writing 'reference/aggMean.html'
# Reading 'man/aggMetadata_parallel_sam.Rd'
# Writing 'reference/aggMetadata_parallel_sam.html'
# Reading 'man/aggMetadata_sam.Rd'
# Writing 'reference/aggMetadata_sam.html'
# Reading 'man/aggregateFeatures_sam.Rd'
# Writing 'reference/aggregateFeatures_sam.html'
# Reading 'man/aggregate_with_matAdj.Rd'
# Writing 'reference/aggregate_with_matAdj.html'
# Reading 'man/aggSum.Rd'
# Writing 'reference/aggSum.html'
# Reading 'man/aggTopn.Rd'
# Writing 'reference/aggTopn.html'
# Reading 'man/barplotEnrichGO_HC.Rd'
# Writing 'reference/barplotEnrichGO_HC.html'
# Reading 'man/barplotGroupGO_HC.Rd'
# Writing 'reference/barplotGroupGO_HC.html'
# Reading 'man/boxPlotD_HC.Rd'
# Writing 'reference/boxPlotD_HC.html'
# Reading 'man/buildGraph.Rd'
# Writing 'reference/buildGraph.html'
# Reading 'man/BuildPalette.Rd'
# Writing 'reference/BuildPalette.html'
# Reading 'man/check.conditions.Rd'
# Writing 'reference/check.conditions.html'
# Reading 'man/checkClusterability.Rd'
# Writing 'reference/checkClusterability.html'
# Reading 'man/CheckDesign.Rd'
# Writing 'reference/CheckDesign.html'
# Reading 'man/classic1wayAnova.Rd'
# Writing 'reference/classic1wayAnova.html'
# Reading 'man/compareNormalizationD_HC.Rd'
# Writing 'reference/compareNormalizationD_HC.html'
# Reading 'man/compute.group.t.test.Rd'
# Writing 'reference/compute.group.t.test.html'
# Reading 'man/compute.t.test.Rd'
# Writing 'reference/compute.t.test.html'
# Reading 'man/ComputeConnexComposants.Rd'
# Writing 'reference/ComputeConnexComposants.html'
# Reading 'man/convertMSnset2QFeatures.Rd'
# Writing 'reference/convertMSnset2QFeatures.html'
# Reading 'man/corrMatrixD_HC.Rd'
# Writing 'reference/corrMatrixD_HC.html'
# Reading 'man/createQFeatures.Rd'
# Writing 'reference/createQFeatures.html'
# Reading 'man/CVDistD_HC.Rd'
# Writing 'reference/CVDistD_HC.html'
# Reading 'man/dapar_hc_chart.Rd'
# Writing 'reference/dapar_hc_chart.html'
# Reading 'man/dapar_hc_ExportMenu.Rd'
# Writing 'reference/dapar_hc_ExportMenu.html'
# Reading 'man/densityPlotD_HC.Rd'
# Writing 'reference/densityPlotD_HC.html'
# Reading 'man/diffAnaComputeFDR.Rd'
# Writing 'reference/diffAnaComputeFDR.html'
# Reading 'man/diffAnalysis.Rd'
# Writing 'reference/diffAnalysis.html'
# Reading 'man/diffAnaVolcanoplot_rCharts.Rd'
# Writing 'reference/diffAnaVolcanoplot_rCharts.html'
# Reading 'man/diff_analysis_sam.Rd'
# Writing 'reference/diff_analysis_sam.html'
# Reading 'man/display.CC.visNet.Rd'
# Writing 'reference/display.CC.visNet.html'
# Reading 'man/enrich_GO.Rd'
# Writing 'reference/enrich_GO.html'
# Reading 'man/filterFeaturesSam.Rd'
# Writing 'reference/filterFeaturesSam.html'
# Reading 'man/find_MEC_matrix.Rd'
# Writing 'reference/find_MEC_matrix.html'
# Reading 'man/formatLimmaResult.Rd'
# Writing 'reference/formatLimmaResult.html'
# Reading 'man/formatPHResults.Rd'
# Writing 'reference/formatPHResults.html'
# Reading 'man/fudge2LRT.Rd'
# Writing 'reference/fudge2LRT.html'
# Reading 'man/get.pep.prot.cc.Rd'
# Writing 'reference/get.pep.prot.cc.html'
# Reading 'man/GetAdjMat.Rd'
# Writing 'reference/GetAdjMat.html'
# Reading 'man/GetDetailedNbPeptides.Rd'
# Writing 'reference/GetDetailedNbPeptides.html'
# Reading 'man/GetDetailedNbPeptidesUsed.Rd'
# Writing 'reference/GetDetailedNbPeptidesUsed.html'
# Reading 'man/getListNbValuesInLines.Rd'
# Writing 'reference/getListNbValuesInLines.html'
# Reading 'man/GetNbPeptidesUsed.Rd'
# Writing 'reference/GetNbPeptidesUsed.html'
# Reading 'man/getQuantile4Imp.Rd'
# Writing 'reference/getQuantile4Imp.html'
# Reading 'man/Get_AllComparisons.Rd'
# Writing 'reference/Get_AllComparisons.html'
# Reading 'man/GlobalQuantileAlignment.Rd'
# Writing 'reference/GlobalQuantileAlignment.html'
# Reading 'man/GOAnalysisSave.Rd'
# Writing 'reference/GOAnalysisSave.html'
# Reading 'man/GraphPepProt_hc.Rd'
# Writing 'reference/GraphPepProt_hc.html'
# Reading 'man/groupttest.Rd'
# Writing 'reference/groupttest.html'
# Reading 'man/group_GO.Rd'
# Writing 'reference/group_GO.html'
# Reading 'man/hc_logFC_DensityPlot.Rd'
# Writing 'reference/hc_logFC_DensityPlot.html'
# Reading 'man/hc_mvTypePlot2.Rd'
# Writing 'reference/hc_mvTypePlot2.html'
# Reading 'man/heatmap.DAPAR.Rd'
# Writing 'reference/heatmap.DAPAR.html'
# Reading 'man/heatmapD.Rd'
# Writing 'reference/heatmapD.html'
# Reading 'man/histPValue_HC.Rd'
# Writing 'reference/histPValue_HC.html'
# Reading 'man/HypothesisTestMethods.Rd'
# Writing 'reference/HypothesisTestMethods.html'
# Reading 'man/imputeMethodsDapar.Rd'
# Writing 'reference/imputeMethodsDapar.html'
# Reading 'man/impute_dapar.Rd'
# Writing 'reference/impute_dapar.html'
# Reading 'man/impute_det_quant.Rd'
# Writing 'reference/impute_det_quant.html'
# Reading 'man/impute_fixed_value.Rd'
# Writing 'reference/impute_fixed_value.html'
# Reading 'man/impute_knn_by_conditions.Rd'
# Writing 'reference/impute_knn_by_conditions.html'
# Reading 'man/impute_matrix_dapar.Rd'
# Writing 'reference/impute_matrix_dapar.html'
# Reading 'man/impute_MEC_Methods.Rd'
# Writing 'reference/impute_MEC_Methods.html'
# Reading 'man/impute_mi.Rd'
# Writing 'reference/impute_mi.html'
# Reading 'man/impute_mle_dapar.Rd'
# Writing 'reference/impute_mle_dapar.html'
# Reading 'man/impute_pa.Rd'
# Writing 'reference/impute_pa.html'
# Reading 'man/impute_pa2.Rd'
# Writing 'reference/impute_pa2.html'
# Reading 'man/impute_POV_Methods.Rd'
# Writing 'reference/impute_POV_Methods.html'
# Reading 'man/impute_slsa.Rd'
# Writing 'reference/impute_slsa.html'
# Reading 'man/inner.aggregate.iter.Rd'
# Writing 'reference/inner.aggregate.iter.html'
# Reading 'man/inner.aggregate.topn.Rd'
# Writing 'reference/inner.aggregate.topn.html'
# Reading 'man/inner.mean.Rd'
# Writing 'reference/inner.mean.html'
# Reading 'man/inner.sum.Rd'
# Writing 'reference/inner.sum.html'
# Reading 'man/is.MV.Rd'
# Writing 'reference/is.MV.html'
# Reading 'man/is.OfType.Rd'
# Writing 'reference/is.OfType.html'
# Reading 'man/LH0.lm.Rd'
# Writing 'reference/LH0.lm.html'
# Reading 'man/LH0.Rd'
# Writing 'reference/LH0.html'
# Reading 'man/LH1.lm.Rd'
# Writing 'reference/LH1.lm.html'
# Reading 'man/LH1.Rd'
# Writing 'reference/LH1.html'
# Reading 'man/limma.complete.test.Rd'
# Writing 'reference/limma.complete.test.html'
# Reading 'man/listSheets.Rd'
# Writing 'reference/listSheets.html'
# Reading 'man/LOESS.Rd'
# Writing 'reference/LOESS.html'
# Reading 'man/make.contrast.Rd'
# Writing 'reference/make.contrast.html'
# Reading 'man/make.design.1.Rd'
# Writing 'reference/make.design.1.html'
# Reading 'man/make.design.2.Rd'
# Writing 'reference/make.design.2.html'
# Reading 'man/make.design.3.Rd'
# Writing 'reference/make.design.3.html'
# Reading 'man/make.design.Rd'
# Writing 'reference/make.design.html'
# Reading 'man/matAdjStats.Rd'
# Writing 'reference/matAdjStats.html'
# Reading 'man/MeanCentering.Rd'
# Writing 'reference/MeanCentering.html'
# Reading 'man/mvHisto_HC.Rd'
# Writing 'reference/mvHisto_HC.html'
# Reading 'man/mvImage.Rd'
# Writing 'reference/mvImage.html'
# Reading 'man/mvPerLinesHistoPerCondition_HC.Rd'
# Writing 'reference/mvPerLinesHistoPerCondition_HC.html'
# Reading 'man/mvPerLinesHisto_HC.Rd'
# Writing 'reference/mvPerLinesHisto_HC.html'
# Reading 'man/MVrowsTagToOne.Rd'
# Writing 'reference/MVrowsTagToOne.html'
# Reading 'man/MVrowsTagToOne_HB.Rd'
# Writing 'reference/MVrowsTagToOne_HB.html'
# Reading 'man/nEmptyLines.Rd'
# Writing 'reference/nEmptyLines.html'
# Reading 'man/nonzero.Rd'
# Writing 'reference/nonzero.html'
# Reading 'man/normalizeD.Rd'
# Writing 'reference/normalizeD.html'
# Reading 'man/normalizeMethods.dapar.Rd'
# Writing 'reference/normalizeMethods.dapar.html'
# Reading 'man/pepa.test.Rd'
# Writing 'reference/pepa.test.html'
# Reading 'man/plotJitterCC.Rd'
# Writing 'reference/plotJitterCC.html'
# Reading 'man/plotJitter_hc.Rd'
# Writing 'reference/plotJitter_hc.html'
# Reading 'man/plotPCA_Eigen_hc.Rd'
# Writing 'reference/plotPCA_Eigen_hc.html'
# Reading 'man/plotPCA_Ind.Rd'
# Writing 'reference/plotPCA_Ind.html'
# Reading 'man/plotPCA_Var.Rd'
# Writing 'reference/plotPCA_Var.html'
# Reading 'man/postHocTest.Rd'
# Writing 'reference/postHocTest.html'
# Reading 'man/POV_impute_det_quant.Rd'
# Writing 'reference/POV_impute_det_quant.html'
# Reading 'man/POV_impute_knn_by_conditions.Rd'
# Writing 'reference/POV_impute_knn_by_conditions.html'
# Reading 'man/POV_impute_slsa.Rd'
# Writing 'reference/POV_impute_slsa.html'
# Reading 'man/proportionConRev_HC.Rd'
# Writing 'reference/proportionConRev_HC.html'
# Reading 'man/QuantileCentering.Rd'
# Writing 'reference/QuantileCentering.html'
# Reading 'man/readExcel.Rd'
# Writing 'reference/readExcel.html'
# Reading 'man/removeAdditionalCol.Rd'
# Writing 'reference/removeAdditionalCol.html'
# Reading 'man/restore_MEC_matrix.Rd'
# Writing 'reference/restore_MEC_matrix.html'
# Reading 'man/rowdata_stats_Aggregation_sam.Rd'
# Writing 'reference/rowdata_stats_Aggregation_sam.html'
# Reading 'man/samLRT.Rd'
# Writing 'reference/samLRT.html'
# Reading 'man/scatterplotEnrichGO_HC.Rd'
# Writing 'reference/scatterplotEnrichGO_HC.html'
# Reading 'man/setMEC.Rd'
# Writing 'reference/setMEC.html'
# Reading 'man/standardise.Rd'
# Writing 'reference/standardise.html'
# Reading 'man/standardiseMeanIntensities.Rd'
# Writing 'reference/standardiseMeanIntensities.html'
# Reading 'man/SumByColumns.Rd'
# Writing 'reference/SumByColumns.html'
# Reading 'man/test.design.Rd'
# Writing 'reference/test.design.html'
# Reading 'man/translatedRandomBeta.Rd'
# Writing 'reference/translatedRandomBeta.html'
# Reading 'man/t_test_sam.Rd'
# Writing 'reference/t_test_sam.html'
# Reading 'man/univ_AnnotDbPkg.Rd'
# Writing 'reference/univ_AnnotDbPkg.html'
# Reading 'man/violinPlotD.Rd'
# Writing 'reference/violinPlotD.html'
# Reading 'man/visualizeClusters.Rd'
# Writing 'reference/visualizeClusters.html'
# Reading 'man/vsn.Rd'
# Writing 'reference/vsn.html'
# Reading 'man/wrapper.pca.Rd'
# Writing 'reference/wrapper.pca.html'
# Reading 'man/wrapperClassic1wayAnova.Rd'
# Writing 'reference/wrapperClassic1wayAnova.html'
# Reading 'man/wrapperRunClustering.Rd'
# Writing 'reference/wrapperRunClustering.html'
# Reading 'man/writeQFeaturesToExcel.Rd'
# Writing 'reference/writeQFeaturesToExcel.html'
# -- Previewing site -------------------------------------------------------------
#   Warning messages:
# 1: multiple methods tables found for ‘rowRanges’ 
# 2: replacing previous import ‘dplyr::last’ by ‘data.table::last’ when loading ‘DAPAR2’ 
# 3: package ‘testthat’ was built under R version 4.0.2


#############################################################################################
# copy paste vignettes/ from master/Prostar for test
build_articles_index()
# Writing 'articles/index.html'
# [1] TRUE
# Warning message:
#   In '_pkgdown.yml', topic must be valid R code
# x Not '``' 

# need _pkgdown.yml. Next to build_tutorials(). build_articles_index() after build_site
build_articles_index() # do nothing


#############################################################################################
build_news()
# -- Building news -----------------------------------------------------------------------------------
#   Writing 'news/index.html'
# -- Previewing site ---------------------------------------------------------------------------------



#################################       END       #################################
build_site()
