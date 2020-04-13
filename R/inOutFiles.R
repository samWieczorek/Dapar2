#' Sets the MEC tag in the origin parameter
#' @title Sets the MEC tag in the OriginOfValues
#' @param origin xxxxx
#' @param qData xxxxxx
#' @param conds xxxxxxx
#' @return An instance of class \code{DataFrame} which is an update of the parameter names origin
#' @author Samuel Wieczorek
#' @examples 
#' @export
setMEC <- function(origin, qData,conds){
  
  #if (is.null( obj@experimentData@other$OriginOfValues)){return()}
  
  u.conds <- unique(conds)
  nbCond <- length(u.conds)
  
  for (c in 1:nbCond){
    ind <- which(conds == u.conds[c])
    if (length(ind) == 1) {
      lst.NA <- which(is.na(qData[,ind]))
    } else {
      lst.NA <- which(apply(is.na(qData[,ind]), 1, sum)==length(ind))
    }
    if (length(lst.NA) > 0)
    {
      origin[lst.NA,ind] <- "MEC"
    }
  }
  return(origin)
}




#' Sets the MEC tag in the origin parameter
#' @title Sets the MEC tag in the OriginOfValues
#' @param obj An object of class \code{Features}
#' @param i The indice of the assay in obj to be used to build the Origin of values. Generally, it is the first one
#' @param namesOrigin A list of names which are known to contain method of identification of peptides.Thoses names must figure in
#' the colnames of the rowData of an assay
#' @return An instance of class \code{DataFrame} Which has the same dimensions as the quantitaive datas (assay) ans contains the origin
#' of values with tags such as MEC, POV
#' @author Samuel Wieczorek
#' @examples 
#' @export
addOriginOfValues <- function(obj, i, namesOrigin=NULL){
  
  if ( ncol(assay(obj,i)) != length(namesOrigin)){
    warning('eee')
    return(NULL)
  }
  
  if (!is.null(namesOrigin))
  {
    OriginOfValues <- rowData(obj[[i]])[,namesOrigin]
  } else {   
    OriginOfValues <- DataFrame(matrix(rep("unknown", nrow(assay(obj,i))*nrow(colData(obj))), 
                                        nrow=nrow(assay(obj,i)),
                                        ncol=nrow(colData(obj))))
  }
  
  ## Identification of all NA as POV
  tmp <- as.matrix(OriginOfValues)
  sel <- is.na(assay(obj,i))
  if (sum(sel) > 0){
    tmp[sel] <- 'POV'
    OriginOfValues <- as(tmp, 'DataFrame')
  }
  
  #rownames(OriginOfValues) <- rownames(Biobase::fData(obj))
  colnames(OriginOfValues) <- paste0("OriginOfValue_",rownames(colData(obj)))
  colnames(OriginOfValues) <- gsub(".", "_", colnames(OriginOfValues), fixed=TRUE)
  
  
  OriginOfValues <- setMEC(origin = OriginOfValues, 
                           qData = assay(obj[['original']]), 
                           conds = colData(obj)$Condition)
  
  
  return(OriginOfValues)
}






#' Builds an object of class \code{Features} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description.
#' 
#' @title Creates an object of class \code{MSnSet} from text file
#' @param data The name of a tab-separated file that contains the data.
#' @param sample A dataframe describing the samples (in lines).
#' @param indExpData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{fData()} table
#' of the \code{MSnSet} object.
#' @param keyId The indice of the column containing the ID of entities 
#' (peptides or proteins)
#' @param namesOrigin xxxxxxxxxxx
#' @param logTransform A boolean value to indicate if the data have to be
#' log-transformed (Default is FALSE)
#' @param forceNA A boolean value to indicate if the 0 and NaN values of
#' intensity have to be replaced by NA (Default is FALSE)
#' @param typeOfData A string that indicates whether the dataset is about
#' @param parentProtId For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' proteins and is used to generate adjacency matrix and process to aggregation.
#' @return An instance of class \code{Features}.
#' @author Samuel Wieczorek
#' @examples 
#' data.file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' data <- read.table(exprsFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' sample <- read.table(sampleFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' indExpData <- c(56:61)
#' namesOrigin <- colnames(data)[43:48]
#' parentId <- 'Protein_group_IDs'
#' keyid <- 'Sequence'
#' ft <- createFeatures(data, sample,indExpData,  keyId = keyid, namesOrigin = namesOrigin, typeOfData = "peptide", parentProtId = parentId, forceNA=TRUE)
#' @import Features
#' @export
createFeatures <- function(data,
                           sample=NULL,
                           indExpData,
                           keyId=NULL,
                           namesOrigin = NULL,
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
  
  # As the function addOribingOfValues is based on the presence of NA in quanti data,
  # if forceNA is not set as TRUE, the previous function cannot ben run
  origin <- DataFrame()
  if (isTRUE(forceNA)) {
    obj <- zeroIsNA(obj,seq_along(obj))
    origin <- addOriginOfValues(obj, 1, namesOrigin)
    metadata(obj)$OriginOfValues <- colnames(origin)
    rowData(obj[['original']]) <- cbind(rowData(obj[['original']]), origin)

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
                        RawPValues = FALSE,
                        OriginOfValues = colnames(origin)
  )
 
  return(obj)
}

