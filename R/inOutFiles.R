

#' Builds an object of class \code{MSnSet} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description. It differs from
#' the original \code{MSnSet} builder which requires three separated files 
#' tabulated-like quantitative proteomic data into a \code{MSnSet} object,
#' including metadata.
#' 
#' @title Creates an object of class \code{MSnSet} from text file
#' @param file The name of a tab-separated file that contains the data.
#' @param metadata A dataframe describing the samples (in lines).
#' @param indExpData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{fData()} table
#' of the \code{MSnSet} object.
#' @param indFData The name of column in \code{file} that will be the name of
#' rows for the \code{exprs()} and \code{fData()} tables
#' @param keyId The indice of the column containing the ID of entities 
#' (peptides or proteins)
#' @param indexForOriginOfValue xxxxxxxxxxx
#' @param logData A boolean value to indicate if the data have to be
#' log-transformed (Default is FALSE)
#' @param replaceZeros A boolean value to indicate if the 0 and NaN values of
#' intensity have to be replaced by NA (Default is FALSE)
#' @param typeOfData A string that indicates whether the dataset is about
#' @param parentProtId xxxx
#' @param versions A list of the following items: Prostar_Version, DAPAR_Version
#' peptides or proteins.
#' @return An instance of class \code{MSnSet}.
#' @author Florence Combes, Samuel Wieczorek
#' @examples 
#' require(Matrix)
#' data.file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' data <- read.table(exprsFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' sample <- read.table(sampleFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' indExpData <- c(56:61)
#' keyid <- 'Sequence'
#' createMSnset(data, sample,indExpData,  keyid, indexForOriginOfValue = NULL, typeOfData = "peptide")
#' @importFrom Features
#' @export
createFeatures <- function(data,
                           sample=NULL,
                           indExpData,
                           keyId=NULL,
                           indexForOriginOfValue = NULL,
                           logTransform=FALSE, 
                           forceNA=FALSE,
                           typeOfData=NULL,
                           parentProtId = NULL){
  
  
  if (is.null(keyId)) {
    obj <- readFeatures(data, 
                        ecol=indExpData, 
                        name='original')
  } else {
    obj <- readFeatures(data, 
                        ecol=indExpData, 
                        name='original',
                        fnames = keyId)
  }
  
  
  ## Encoding the sample data
  sample <- lapply(sample,function(x){ gsub(".", "_", x, fixed=TRUE)})
  colData(obj)@listData <- sample
  
  
  ## Replace all '.' by '_' in names
  #colnames(fd) <- gsub(".", "_", colnames(data)[indFData], fixed=TRUE)
  #colnames(Intensity) <- gsub(".", "_", colnames(data)[indExpData], fixed=TRUE)
  
  
  if (isTRUE(forceNA)) {
    obj <- zeroIsNA(obj,seq_along(obj))
  }
  
  
  if (isTRUE(logTransform)) {
    obj <- addAssay(obj, logTransform(obj[['original']]),name = "original_log")
    obj <- addAssayLinkOneToOne(obj, from = "original", to = "original_log")
  }
  
  
  daparVersion <- if (length(grep('DAPAR', installed.packages())) > 0) {
    installed.packages()["DAPAR","Version"]
  } else {
    'NA'
  }
  
  ProstarVersion <- if (length(grep('Prostar2', installed.packages())) > 0) {
    installed.packages()["Prostar2","Version"]
  } else {
    'NA'
  }
  
  metadata(obj) <- list(versions = list(Prostar_Version = ProstarVersion,
                                        DAPAR_Version = daparVersion),
                        parentProtId = parentProtId,
                        params = list(),
                        typeOfData = typeOfData,
                        originOfValues = NULL,
                        RawPValues = FALSE
  )
  
  
  obj <- addOriginOfValue(obj,indexForOriginOfValue)
  
  return(obj)
}

