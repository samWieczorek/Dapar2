

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
#' write2excel(ft, 'foo.xlsx')
#' write2excel(ft, 1, 'foo.xlsx')
#' 
#' #---------------------------------------
#' # Export one assay
#' #---------------------------------------
#' 
#' 
#' obj <- Exp1_R25_pept[seq_len(1000)]
#' writeQFeaturesToExcel(obj, "foo")
#' }
#' 
#' @export
#' 
#' 
#' @rdname QFeatures-excel
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
                          exp.design(object),
                          filename, 
                          ...)
              
              
            }
            
  })


##' @exportMethod write2excel
##' @rdname QFeatures-excel
setMethod("write2excel", "SummarizedExperiment",
          function(object, ...)
            .write2excel(object, exp.design, filename, ...))







#' @importFrom openxlsx createStyle createWorkbook addWorksheet writeData addStyle writeData
.write2excel <- function(object, exp.design, filename) {
  name <- paste(filename, ".xlsx", sep="")
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
    
    
    unique.tags <- NULL
    if (!is.null(tags) && !is.null(colors)){
      unique.tags <- unique(as.vector(as.matrix(tags)))
      cond.colors <- sum(unique.tags %in% names(colors)) == length(unique.tags)
      
      if (cond.colors){
        lapply(1:length(colors), function(x){
          list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
          openxlsx::addStyle(wb,
                             sheet = i.sheet,
                             cols = list.tags[ ,"col"],
                             rows = list.tags[ ,"row"] + 1, 
                             style = openxlsx::createStyle(fgFill = colors[x])
                             )
          })
      } else {
        warning("The length of colors vector must be equal to the 
        number of different tags. As is it not the case, colors are ignored")
      }
    }
    
    
    # Write the rowData table to the second sheet
    i.sheet <- 2
    
    # Write only one-dimensional slots
    
    oneDim <- which(lapply(rowData(object), is.vector) == 1)
    openxlsx::addWorksheet(wb, "Samples Meta Data")
    openxlsx::writeData(wb, 
                        sheet = i.sheet, 
                        rowData(object)[, oneDim], 
                        rowNames = FALSE)
    
    
    # Add colors for sample data sheet
    u_conds <- unique(Biobase::pData(obj)$Condition)
    colors <- setNames(DAPAR::ExtendPalette(length(u_conds)),
                       u_conds)
    colors[['blank']] <- 'white'
    
    tags <- Biobase::pData(obj)
    tags[,] <- 'blank'
    tags$Sample.name <- Biobase::pData(obj)$Condition
    tags$Condition <- Biobase::pData(obj)$Condition
    
    unique.tags <- NULL
    if (!is.null(tags) && !is.null(colors)){
      unique.tags <- unique(as.vector(as.matrix(tags)))
      if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
        warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
      if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
        lapply(1:length(colors), function(x){
          list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
          openxlsx::addStyle(wb,
                             sheet = n,
                             cols = list.tags[ ,"col"],
                             rows = list.tags[ ,"row"] + 1, 
                             style = openxlsx::createStyle(fgFill = colors[x])
          )
        })
      }
    }
    
    
    ## Add the experimental design to the third sheet
    
    n <- 3
    if (dim(Biobase::fData(obj))[2] != 0){
      openxlsx::addWorksheet(wb, "Feature Meta Data")
      openxlsx::writeData(wb, 
                          sheet = n, 
                          cbind(ID = rownames(Biobase::fData(obj)),
                                Biobase::fData(obj)), rowNames = FALSE)
    }
    
    colors <- as.list(setNames(mc$color, mc$node))
    tags <- cbind(keyId = rep('identified', nrow(obj)),
                  Biobase::fData(obj)
    )
    
    tags[,] <- 'identified'
    tags[, 1 + which(colnames(fData(obj)) %in% obj@experimentData@other$names_metacell)] <- GetMetacell(obj)
    
    unique.tags <- NULL
    if (!is.null(tags) && !is.null(colors)){
      unique.tags <- unique(as.vector(as.matrix(tags)))
      if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
        warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
      if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
        lapply(1:length(colors), function(x){
          list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
          openxlsx::addStyle(wb,
                             sheet = n,
                             cols = list.tags[ ,"col"],
                             rows = list.tags[ ,"row"] + 1, 
                             style = openxlsx::createStyle(fgFill = colors[x])
          )
        })
      }
    }
    
    # Add GO tab
    if (!is.null(obj@experimentData@other$GGO_analysis))
    {
      l <- length(obj@experimentData@other$GGO_analysis$ggo_res)
      for (i in 1:l){
        n <- n +1
        level <- as.numeric(obj@experimentData@other$GGO_analysis$levels[i])
        openxlsx::addWorksheet(wb, paste("Group GO - level ", level, sep=""))
        openxlsx::writeData(wb, sheet=n, obj@experimentData@other$GGO_analysis$ggo_res[[i]]$ggo_res@result)
      }
    }
    
    
    
    if (!is.null(obj@experimentData@other$EGO_analysis))
    {
      n <- n +1
      openxlsx::addWorksheet(wb, "Enrichment GO")
      openxlsx::writeData(wb, sheet=n, obj@experimentData@other$EGO_analysis$ego_res@result)
      
    }
    
    openxlsx::saveWorkbook(wb, name, overwrite=TRUE)
    return(name)
    
    
  }
  
}

