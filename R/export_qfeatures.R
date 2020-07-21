

#' @title This function exports a \code{QFeatures} object to a Excel file.
#' 
#' @description This function exports a \code{QFeatures} data object to a Excel file.
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
#' Sys.setenv("R_ZIPCMD"= Sys.which("zip"))
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept[1:1000]
#' writeQFeaturesToExcel(obj, "foo")
#' }
#' 
#' @export
#' 
#' @importFrom openxlsx createStyle createWorkbook addWorksheet writeData addStyle writeData
#' 
writeQFeaturesToExcel <- function(obj, filename) {
  
  POV_Style <- openxlsx::createStyle(fgFill = "lightblue")
  MEC_Style <- openxlsx::createStyle(fgFill = "orange")
  
  name <- paste(filename, ".xlsx", sep="")
  wb <- openxlsx::createWorkbook(name)
  
  
  #
  # Add quantitative tabs for each assay
  #
  for (i in 1:length(obj)){
    openxlsx::addWorksheet(wb, paste0("Quanti Data for ", names(obj)[i]))
    openxlsx::writeData(wb, 
                        sheet=i,
                        cbind(ID = rownames(assay(obj, i)),
                              assay(obj, i)),
                        rowNames = FALSE)
    
    
    if (is.null(metadata(obj)$OriginOfValues)){
      listPOV <-  which(is.na(assay(obj, i)), arr.ind=TRUE)
    } else {
      mat <- as.data.frame(rowData(obj[[i]])[,metadata(obj)$OriginOfValues])
      listPOV <- which(mat=="POV", arr.ind=TRUE)
      listMEC <- which(mat=="MEC", arr.ind=TRUE)
    }
    
    openxlsx::addStyle(wb, sheet=i, cols = listPOV[,"col"]+1, rows = listPOV[,"row"]+1, style = POV_Style)
    openxlsx::addStyle(wb, sheet=i, cols = listMEC[,"col"]+1, rows = listMEC[,"row"]+1, style = MEC_Style)
  }
  
  
  #
  # Add sample tab
  #
  openxlsx::addWorksheet(wb, "Samples Meta Data")
  openxlsx::writeData(wb, sheet=(1 + length(obj)), colData(obj), rowNames = FALSE)
  
  
  #
  # Add quantitative tabs for each assay
  #
  for (i in 1:length(obj)){
    offset <- 1 + length(obj) + i
    openxlsx::addWorksheet(wb, paste0("Meta Data for ", names(obj)[i]))
    
    openxlsx::writeData(wb, sheet=offset, cbind(ID = rownames(assay(obj, i)),
                                                rowData(obj[[i]])), rowNames = FALSE)
  }
  
  # if (!is.null(obj@experimentData@other$GGO_analysis))
  # {
  #   l <- length(obj@experimentData@other$GGO_analysis$ggo_res)
  #   for (i in 1:l){
  #     n <- n +1
  #     level <- as.numeric(obj@experimentData@other$GGO_analysis$levels[i])
  #     openxlsx::addWorksheet(wb, paste("Group GO - level ", level, sep=""))
  #     openxlsx::writeData(wb, sheet=n, obj@experimentData@other$GGO_analysis$ggo_res[[i]]$ggo_res@result)
  #   }
  # }
  # 
  # if (!is.null(obj@experimentData@other$EGO_analysis))
  # {
  #   n <- n +1
  #   openxlsx::addWorksheet(wb, "Enrichment GO")
  #   openxlsx::writeData(wb, sheet=n, obj@experimentData@other$EGO_analysis$ego_res@result)
  #   
  # }
  
  openxlsx::saveWorkbook(wb, name, overwrite=TRUE)
  return(name)
  
}


#' @title This function reads a sheet of an Excel file and put the data into a data.frame.
#' 
#' @param file The name of the Excel file.
#' 
#' @param extension The extension of the file
#' 
#' @param sheet The name of the sheet
#' 
#' @return A data.frame
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @importFrom readxl read_excel
#' 
readExcel <- function(file, extension, sheet){
  # data <- NULL
  # if (extension=="xls") {
  #     data <- readxl::read_xls(file, sheet)
  # }
  # else if (extension=="xlsx") {
  #     data <- readxl::read_xlsx(file, sheet)
  # }
  # return(as.data.frame(data,asIs=T))
  
  #options(digits=10)
  data <- NULL
  data <- readxl::read_excel(file, sheet)
  
  return(as.data.frame(data,asIs=T, stringsAsFactors=F))
  
}



#' @title This function returns the list of the sheets names in a Excel file.
#' 
#' @param file The name of the Excel file.
#' 
#' @return A vector
#' 
#' @author Samuel Wieczorek
#'  
#' @export
#' 
#' @importFrom openxlsx getSheetNames
#' 
listSheets <- function(file){
  ll <- openxlsx::getSheetNames(file)
  return(ll)
  
}



#' #' @title Exports a MSnset dataset into a zip archive containing three zipped CSV files.
#' #' 
#' #' @param obj An object of class \code{Features}.
#' #' 
#' #' @param fname The name of the archive file.
#' #' 
#' #' @return A compressed file
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples
#' #' \donttest{
#' #' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' #' obj <- Exp1_R25_pept[1:1000]
#' #' writeMSnsetToCSV(obj, "foo")
#' #' }
#' #' 
#' #' @export
#' #' 
#' #' @import zip
#' #' 
#' writeMSnsetToCSV <- function(obj, fname){
#'   
#'   write.csv(Biobase::exprs(obj), paste(tempdir(), "exprs.csv", sep='/'))
#'   write.csv(Biobase::fData(obj), paste(tempdir(), "fData.csv", sep='/'))
#'   write.csv(Biobase::pData(obj), paste(tempdir(), "pData.csv", sep='/'))
#'   files <- c(paste(tempdir(), "exprs.csv", sep='/'),
#'              paste(tempdir(), "fData.csv", sep='/'),
#'              paste(tempdir(), "pData.csv", sep='/'))
#'   zip(fname, files, zip = Sys.getenv("R_ZIPCMD", "zip"))
#'   
#'   return(fname)
#' }