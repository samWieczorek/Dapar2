#' @title This function exports a \code{QFeatures} object to a Excel file.
#' 
#' @description 
#' 
#' This function exports a \code{QFeatures} data object to a Excel file.
#' Each of the three data.frames in the \code{QFeatures} object (ie expression data,
#' colData and rowData are respectively integrated into separate sheets in
#' the Excel file). The colored cells in the experimental data correspond to the 
#' original missing values which have been imputed.
#' 
#' @param obj An object of class \code{QFeatures}.
#' 
#' @param filename A character string for the name of the Excel file.
#' 
#' @return A Excel file (.xlsx)
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
#' 
#' @export
#' 
#' @name QFeatures-excel
#'
NULL


##' @exportMethod write2excel
##' @rdname QFeatures-excel
setMethod("write2excel", "QFeatures",
          function(object, 
                   i = NULL,
                   filename = "newFile", ...) {
            if (isEmpty(object))
              return(NULL)
            
            
            if (is.null(i)){
              # One exports all the QFeatures object
              
            } else {
              # One exports only one SE
              write2excel(object[[i]],
                          filename, 
                          design(object),
                          ...)
              
              
            }
            
          })


##' @exportMethod write2excel
##' @rdname QFeatures-excel
setMethod("write2excel", "SummarizedExperiment",
          function(object, filename, exp.design, ...)
            .write2excel(object, filename, exp.design, ...))







#' @importFrom openxlsx createStyle createWorkbook addWorksheet writeData addStyle writeData
.write2excel <- function(object, exp.design, fname) {
  name <- paste0(fname, ".xlsx", sep="")
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


#' @param wb xxx
#' @param n xxx
#' @param tags xxx
#' @param colors xxx
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