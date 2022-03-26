#' @title Exports a \code{QFeatures} object to a Excel file.
#' 
#' @description 
#' 
#' This function exports an instance of the class `QFeatures` to a Excel file.
#' The resulting file is composed of four sheets:
#' 
#' - `quantitative data` which contains the content of [assay()] object whith a 
#' color code for each cell w.r.t. to cell quantitative metadata.
#' 
#' - `metadata` which is the content of [rowData()] with only one-dimensionnal
#' data (i.e. the adjacencyMatrix and the qMetadata slots are not part of
#' the sheet),
#' 
#' - `exp. design` which is the content of [colData()]. Each condition in the table
#' is colored with a different color,
#' 
#' - `quantitative metadata` which is the content of [qMetadata()]. There is a color
#' code for the different tags.
#' 
#' 
#' @param object An object of class \code{QFeatures}.
#' @param i xxx
#' @param filename A character string for the name of the Excel file.
#' @param exp.design xxx
#' 
#' @return A Excel file.
#' 
#' @name QFeatures-excel
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \donttest{
#' data(ft)
#' 
#' #---------------------------------------
#' # Export the whole dataset
#' #---------------------------------------
#' 
#' write2excel(ft, filename = 'foo')
#' write2excel(ft, 1, 'foo')
#' }
#'
NULL


#' @exportMethod write2excel
#' @rdname QFeatures-excel
setMethod("write2excel", "QFeatures",
          function(object, 
                   i = NULL,
                   filename = "newFile", ...) {
            if (length(object)==0)
              return(NULL)
            
            
            if (is.null(i)){
              # One exports all the QFeatures object
              
            } else {
              # One exports only one SE
              write2excel(object[[i]],
                          filename, 
                          design.qf(object),
                          ...)
              
              
            }
            
          })


#' @exportMethod write2excel
#' @rdname QFeatures-excel
setMethod("write2excel", "SummarizedExperiment",
          function(object, filename, exp.design, ...)
            .write2excel(object, filename, exp.design, ...)
          )






#' @param object xxx
#' @param filename xxx
#' @param exp.design xxx
#' 
#' @noRd
#' @importFrom stats setNames
.write2excel <- function(object, filename, exp.design) {
  
  if (! requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install openxlsx: BiocManager::install('openxlsx')")
  }
  
  name <- paste0(filename, ".xlsx", sep="")
  wb <- openxlsx::createWorkbook(name)
  # Write the assay data to the first sheet
  i.sheet <- 1
  openxlsx::addWorksheet(wb, "Quantitative data")
  openxlsx::writeData(wb, 
                      sheet = i.sheet, 
                      cbind(ID = rowData(object)[, idcol(object)],
                            assay(object)), 
                      rowNames = FALSE)
  
  
  # Add colors to quantitative table
  mc <- qMetadata.def(typeDataset(object))
  colors <- as.list(setNames(mc$color, mc$node))
  tags <- cbind(keyId = rep('identified', nrow(object)),
                qMetadata(object)
                )
  
  addColors(wb, i.sheet, tags, colors)
  
  # Write the rowData table to the second sheet
  i.sheet <- 2

  # Write only one-dimensional slots
#browser()
  openxlsx::addWorksheet(wb, "Exp. design")
  openxlsx::writeData(wb,
                      sheet = i.sheet,
                      exp.design,
                      rowNames = FALSE)


  # Add colors for sample data sheet
  u_conds <- unique(exp.design$Condition)
  colors <- setNames(ExtendPalette(length(u_conds)),u_conds)
  colors[['blank']] <- 'white'

  tags <- as.data.frame(exp.design)
  tags[,] <- 'blank'
  tags$Sample.name <- exp.design$Condition
  tags$Condition <- exp.design$Condition

  addColors(wb, i.sheet, tags, colors)

  # 
  # ## Add the experimental design to the third sheet

  n <- 3
  oneDim <- which(lapply(rowData(object), is.vector) == 1)
  new.rowData <- rowData(object)[, oneDim]
  openxlsx::addWorksheet(wb, "rowData")
    openxlsx::writeData(wb,
                        sheet = n,
                        cbind(ID = new.rowData[, idcol(object)],
                        as.data.frame(new.rowData)),
                        rowNames = FALSE)

    
    # Add the qMetadata information
    n <- 4
    new.rowData <- qMetadata(object)
    openxlsx::addWorksheet(wb, "qMetadata")
    openxlsx::writeData(wb,
                        sheet = n,
                        cbind(ID = rowData(object)[, idcol(object)],
                        new.rowData),
                        rowNames = FALSE)

    
   colors <- as.list(setNames(mc$color, mc$node))
  tags <- cbind(keyId = rep('identified', nrow(new.rowData)),
                new.rowData
                )

  tags[,] <- 'identified'
  tags[, 1 + seq_len(ncol(new.rowData))] <- new.rowData

  addColors(wb, n, tags, colors)
  


  
  openxlsx::saveWorkbook(wb, name, overwrite=TRUE)
  return(name)
  
  
}


#' @param wb A workbook
#' @param n A `integer(1)` which is the number of sheet in the workbook.
#' @param tags xxx
#' @param colors A `character()` which contains the HEX code for colors. The size
#' of this vector must be the same as the number of tags.
#' 
#' @noRd
#' 
addColors <- function(wb, n, tags, colors){
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)){
    unique.tags <- unique(as.vector(as.matrix(tags)))
    conds.colors <- sum(unique.tags %in% names(colors)) == length(unique.tags)
    
    if (conds.colors){
      lapply(seq_len(length(colors)), function(x){
        list.tags <- which(names(colors)[x] == tags, arr.ind=TRUE)
        openxlsx::addStyle(wb,
                           sheet = n,
                           cols = list.tags[ ,"col"],
                           rows = list.tags[ ,"row"] + 1,
                           style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    } else {
      warning("The length of colors vector must be equal to the number of different tags.
              As is it not the case, colors are ignored")
    }
  }
}