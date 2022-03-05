#' @param object An instance of the class `SummarizedExperiment`
#' 
#' @param exp.design xxx
#' 
#' @param subset A vector of index indicating rows to highlight
#' 
#' @param pal.name xxx
#' 
#' @return A violinplot
#' 
#' 
#' @importFrom vioplot vioplot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot.window axis mtext legend points segments plot.new
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' 
#' @rdname descriptive-statistics
#' 
violinPlotD <- function(object, 
                        exp.design, 
                        subset = NULL, 
                        pal.name){

  
  stopifnot(inherits(object, 'SummarizedExperiment'))
  
  qData <- assay(object)
  legend <- exp.design$Sample.name
  
  if(missing(exp.design))
    stop("'exp.design' is missing.")
  
  # In order to view subset, one need the 'idcol' info
  if (!is.null(subset))
    stopifnot(!is.null(idcol(object)))

  
  myColors <- SampleColors(exp.design$Condition, pal.name)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(qData)+1),
                        ylim = c(min(na.omit(qData)),
                                 max(na.omit(qData)))
                        )
  
  title( ylab="Log (intensity)")
  for (i in seq_len(ncol(qData)))
    vioplot::vioplot(na.omit(qData[ ,i]), 
                     col = myColors[i], 
                     add = TRUE, 
                     at = i)

  graphics::axis(2, 
                 yaxp = c(floor(min(na.omit(qData))), 
                          floor(max(na.omit(qData))), 5),
                 las=1)
  
   graphics::axis(side = 1,
                  at = seq_len(ncol(qData)),
                  labels = legend
                  )

  # Display of rows to highlight (index of row in subset) 
  if(!is.null(subset)){
    idVector <- rowData(object)[ ,idcol(object)]
    pal.tracker <- ExtendPalette(length(subset), "Dark2")
    n=0
    for (i in subset) {
      n = n + 1
      for (c in seq_len(ncol(qData)-1)) {
        graphics::segments(y0 = qData[i, c],
                           y1 = qData[i, c + 1],
                           x0 = c,
                           x1 = c + 1,
                           pch = 16,
                           col = pal.tracker[n],
                           lwd = 2
                           )
        graphics::points(y = qData[i,c],
                         x = c,
                         pch = 16,
                         col = pal.tracker[n]
                         )
      }
      graphics::points(y = qData[i, ncol(qData)],
                       x = ncol(qData),
                       pch = 16,
                       col = pal.tracker[n]
                       )
    }
    graphics::legend("topleft",
                     legend = idVector[subset],
                     lty = 1,
                     lwd = 2,
                     col = pal.tracker,
                     pch = 16,
                     bg = "transparent",
                     bty = "n"
                     )
  }
  
}
