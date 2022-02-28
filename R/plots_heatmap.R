#' @title Heatmap of the quantitative proteomics data
#' 
#' @description 
#' 
#' This function is a wrapper to 'heatmap.2' of the package `gplots`
#' that displays a numeric matrix
#' 
#' @section generic plot:
#' 
#' The function [heatmapD()] 
#' 
#' @section plot for missingvalues:
#' 
#' The function [] is inspired from the function 'heatmap.2' 
#' that displays a numeric matrix. For more information, please refer to the help 
#' of the heatmap.2 function.
#' 
#' @return A heatmap
#' 
#' @author Samuel Wieczorek
#' 
#' @name heatmap_plot
#' 
#' @examples
#' 
#' library(QFeatures)
#' data(ft)
#' conds <- design(ft)$Condition
#'
#' #----------------------------------------------
#' # Plots a heatmap for generic quantitative data
#' #----------------------------------------------
#' 
#' heatmapD(assay(ft, 1), conds)
#' 
#' #----------------------------------------------
#' # Plots a heatmap for missing values visualization
#' #----------------------------------------------
#' 
#' library(QFeatures)
#' data(ft)
#' mv.heatmap(assay(ft, 1))
#' 
#' 
#' #----------------------------------------------
#' # Launches a shiny module
#' #----------------------------------------------
#' 
#' if(interactive()){
#'  library(QFeatures)
#'  library(shiny)
#'  library(DaparToolshed)
#'  data(ft)
#'  ui <- mod_heatmap_plot_ui('plot')
#' 
#'  server <- function(input, output, session) {
#'   conds <- design(ft)$Condition
#'  
#'   mod_heatmap_plot_server('plot',
#'                       obj = reactive({ft[[1]]}),
#'                       conds = reactive({conds}),
#'                       pal.name = reactive({'Dark2'})
#'                       )
#'   }
#'  
#'  shinyApp(ui=ui, server=server)
#' }
NULL






#' @param qData A `matrix` or `array` quantitative values
#' 
#' @param conds A `character()` of the conditions (one condition per sample).
#' 
#' @param distfun A `character(1)` defining the distance used by the clustering 
#' algorithm to compute. Default value is 'euclidean'. See \code{help(heatmap.2)}
#' 
#' @param hclustfun the clustering algorithm used to build the dendrogram.
#' Default value is 'complete'. See \code{help(heatmap.2)}
#' 
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom gplots heatmap.2
#' 
#' @export
#' 
#' @rdname heatmap_plot
#' 
heatmapD <- function(qData, 
                     conds, 
                     distfun = "euclidean", 
                     hclustfun = "complete"){
  #Check parameters
  if (!(distfun %in% c("euclidean", "manhattan"))){
      stop("'distfun' is not correct.")
      return (NULL)
  }

  if (!(hclustfun %in%  c("complete", "ward.D", "average"))){
      stop("'hclustfun' is not correct.")
      return (NULL)
  }
  
  if (length(conds) != ncol(qData))
    stop("The length of 'conds' must be equal to the number of samples.")
  
  # .dendro = "none"
  # if (dendro)
  #   .dendro = "row"
  
  # if (isTRUE(dendro) && getNumberOfEmptyLines(qData) != 0)  {
  #     stop("Your dataset contains empty lines: the dendrogram cannot
  # be computed.
  #         Please filter or impute missing values before.")
  #     return (NULL)
  # }
  # else {
  .data <- matrix(qData,
                  ncol = ncol(qData),
                  byrow = FALSE,
                  dimnames = list(rownames(qData), 
                                  colnames(qData)
                                  )
                  )
  
  #
  #heatmap.color <- grDevices::colorRampPalette(c("green", "red"))(n = 1000)
  
  # samples label color
  x = t(.data)
  x[is.na(x)] <- -1e5
  # dist = dist(x, method=distance)
  # dend1 = hclust(dist, method=cluster) %>% as.dendrogram()
  # 
  # cols_branches <- ExtendPalette(n = length(unique((conds))), base = NULL)
  # 
  # dend1 <- dend1 %>% dendextend::color_branches(
  #   k = length(unique(conds)),
  #   col = cols_branches)
  # 
  # dendextend::labels_colors(dend1) <- dendextend::get_leaves_branches_col(dend1)
  # col_labels <- dendextend::get_leaves_branches_col(dend1)
  heatmap(x, Rowv = FALSE, Colv = NULL, na.rm=FALSE)
  
  #   p <- gplots::heatmap.2(
  #   x = t(.data),
  #   distfun = function(x) {
  #     x[is.na(x)] <- -1e5
  #     dist(x, method=distance)
  #   },
  #   hclustfun = function(x) {
  #     x[is.na(x)] <- -1e5
  #     hclust(x, method=cluster)
  #   },
  #   dendrogram = 'none',
  #   Rowv = FALSE,
  #   Colv = FALSE,
  #   col = grDevices::colorRampPalette(c("green", "red"))(n = 1000) ,
  #   density.info = 'none',
  #   key = TRUE,
  #   trace = "none",
  #   scale = "none",
  #   #srtCol=45,
  #   labCol = "",
  #   margins=c(4,12),
  #   #cexRow=1.5,
  #   cexRow = 1.5 + ncol(.data)*-0.011 ,
  #   keysize = 1.5,
  #   lhei = c(1.5, 9),
  #   lwid = c(1.5, 4),
  #   lmat = rbind(4:3, 2:1),
  #   colRow = labels_colors(dend1)
  # )
  
}



#' @param x A `matrix` or `array` containing the quantitative data.
#' 
#' @param col Colors used for the image. Defaults to heat colors (heat.colors).
#' 
#' @param srtCol Angle of column conds, in degrees from horizontal 
#' 
#' @param labCol Character vectors with column conds to use.
#' 
#' @param labRow Character vectors with row conds to use.
#' 
#' @param key Logical indicating whether a color-key should be shown.
#' 
#' @param key.title Main title of the color key. If set to NA no title will 
#' be plotted.
#' 
#' @param main Main title; default to none.
#' 
#' @param ylab y-axis title; default to none.
#' 
#' @export
#' 
#' @importFrom grDevices heat.colors
#' @importFrom graphics image strwidth strheight axis mtext text title layout par plot.new
#' 
#' @rdname heatmap_plot
#' 
mv.heatmap <- function (x,
                        col = grDevices::heat.colors(100),
                        srtCol=NULL,
                        labCol = NULL,
                        labRow = NULL,
                        key = TRUE, 
                        key.title = NULL,
                        main = NULL,  
                        ylab = NULL) {
  
  stopifnot(inherits(x, 'matrix'))
  
  scale01 <- function(x, low = min(x), high = max(x))
    (x - low)/(high - low)

  
  offsetCol = 0.5
  offsetRow = 0.5
  srtRow <- colRow <- colCol <- xlab <- breaks <- na.rm <- NULL
  key.par = list()
  margins = c(5, 5)
  sepcolor <- na.color <- "white"
  keysize = 1.5
  
  
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  x <- x[seq_len(nr),]
  cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  cexCol = 0.2 + 1/log10(nc)
  cexRow = 0.2 + 1/log10(nr)
  iy <- seq_len(nr)
  breaks <- length(col) + 1
  breaks <- seq(min(x, na.rm = na.rm), 
                max(x, na.rm = na.rm), 
                length = breaks)
  
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  lhei <- c(keysize, 4)
  lwid <- c(keysize, 4)
  lmat <- rbind(4:3, 2:1)
  lmat[is.na(lmat)] <- 0
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  graphics::par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  
  
  graphics::image(seq_len(nc),
                  seq_len(nr), 
                  x, 
                  xlim = 0.5 + c(0, nc), 
                  ylim = 0.5 + c(0, nr),
                  axes = FALSE,
                  xlab = "", 
                  ylab = "",
                  col = col, 
                  breaks = breaks
                  )
  
  
  if (!is.null(labCol)) {
    graphics::axis(1, 
                   seq_len(nc), 
                   labels = labCol, 
                   las = 2, 
                   line = -0.5 + offsetCol,
                   tick = 0, 
                   cex.axis = cexCol, 
                   hadj = NA, 
                   padj = 0)
  }
  else {
    adjCol = c(1, NA)
    xpd.orig <- par("xpd")
    par(xpd = NA)
    xpos <- graphics::axis(1,
                           seq_len(nc),
                           labels = rep("", nc),
                           las = 2, 
                           tick = 0
                           )
    graphics::text(x = xpos, 
                   y = par("usr")[3] - (1 + offsetCol) * strheight("M"), 
                   label = labCol, 
                   adj = adjCol, 
                   cex = cexCol, 
                   srt = srtCol, 
                   col = colCol)
    graphics::par(xpd = xpd.orig)
  }
  
  
  if (!is.null(labRow) ) {
    graphics::axis(4, 
                   iy, 
                   labels = labRow, 
                   las = 5, 
                   line = -0.5 + offsetRow, 
                   tick = 0, 
                   cex.axis = cexRow, 
                   hadj = 0, 
                   padj = NA)
  }
  else {
    xpd.orig <- par("xpd")
    par(xpd = NA)
    ypos <- axis(4, 
                 iy, 
                 labels = rep("", nr), 
                 las = 2, 
                 line = -0.5, 
                 tick = 0
                 )
    graphics::text(x = par("usr")[2] + (1 + offsetRow) * graphics::strwidth("M"), 
                   y = ypos, 
                   labels = labRow, 
                   adj = c(0,NA), 
                   cex = cexRow, 
                   srt = srtRow, 
                   col = colRow
                   )
    graphics::par(xpd = xpd.orig)
  }
  
  
  graphics::par(mar = c(margins[1], 0, 0, 0))
  graphics::plot.new()
  graphics::par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  
  graphics::plot.new()
  if (!is.null(main)) 
    graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
  
  
  if (key) {
    mar <- c(5, 4, 2, 1)
    graphics::par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0) 
      do.call(par, key.par)
    
    tmpbreaks <- breaks
    min.raw <- min.breaks
    max.raw <- max.breaks
    
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
    graphics::image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
                    xaxt = "n", yaxt = "n")
    graphics::par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    xargs <- list(at = xv, labels = lv)
    
    xargs$side <- 1
    do.call(axis, xargs)
    key.xlab <- "Intensity value"
    
    graphics::mtext(side = 1, 
                    key.xlab, 
                    line = par("mgp")[1], 
                    padj = 0.5, 
                    cex = par("cex") * par("cex.lab")
                    )
    
    if (is.null(key.title)) 
      title("Color Key")
  }
  
}


