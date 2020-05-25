#' Sets the MEC tag in the origin parameter
#' 
#' @title Sets the MEC tag in the OriginOfValues
#' 
#' @param origin xxxxx
#' 
#' @param qData xxxxxx
#' 
#' @param conds xxxxxxx
#' 
#' @return An instance of class \code{DataFrame} which is an update of the parameter names origin
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' obj <- Exp1_R25_pept
#' }
#' 
#' @export
#' 
setMEC <- function(origin, qData, conds){
  
  if (missing(origin))
    stop("'origin' is required.")
  
  if (missing(qData))
    stop("'qData' is required.")
  if (missing(conds))
    stop("'conds' is required.")
  

  
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
#' 
#' @title Sets the MEC tag in the OriginOfValues
#' 
#' @param obj An object of class \code{Features}
#' 
#' @param i The indice of the assay in obj to be used to build the Origin of values. Generally, it is the first one
#' 
#' @param namesOrigin A list of names which are known to contain method of identification of peptides.Thoses names must figure in
#' the colnames of the rowData of an assay
#' 
#' @return An instance of class \code{DataFrame} Which has the same dimensions as the quantitaive datas (assay) ans contains the origin
#' of values with tags such as MEC, POV
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata2')
#' origin <- c("Identification_type_C_R1", "Identification_type_C_R2", "Identification_type_C_R3", "Identification_type_D_R1", "Identification_type_D_R2", "Identification_type_D_R3")
#' obj <- addOriginOfValues(Exp1_R25_pept, 2, origin)
#' }
#' 
#' @importFrom methods as
#' 
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' 
#' @export
addOriginOfValues <- function(obj, i, namesOrigin=NULL){
  
  if ( ncol(assay(obj,i)) != length(namesOrigin)){
    stop("The number of samples in the assay must be equal to the length of 'namesOrigin'")
  }
  
  if (!is.null(namesOrigin))
  {
    OriginOfValues <- rowData(obj[[i]])[,namesOrigin]
  } else {   
    OriginOfValues <- MultiAssayExperiment::DataFrame(matrix(rep("unknown", nrow(assay(obj,i))*nrow(SummarizedExperiment::colData(obj))), 
                                        nrow=nrow(assay(obj,i)),
                                        ncol=nrow(SummarizedExperiment::colData(obj))))
  }
  
  ## Identification of all NA as POV
  tmp <- as.matrix(OriginOfValues)
  sel <- is.na(assay(obj,i))
  if (sum(sel) > 0){
    tmp[sel] <- 'POV'
    OriginOfValues <- as(tmp, 'DataFrame')
  }
  
  colnames(OriginOfValues) <- paste0("OriginOfValue_",rownames(SummarizedExperiment::colData(obj)))
  colnames(OriginOfValues) <- gsub(".", "_", colnames(OriginOfValues), fixed=TRUE)
  
  
  OriginOfValues <- setMEC(origin = OriginOfValues, 
                           qData = assay(obj[['original']]), 
                           conds = SummarizedExperiment::colData(obj)$Condition)
  
  
  return(OriginOfValues)
}






#' Builds an object of class \code{Features} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description.
#' 
#' @title Creates an object of class \code{MSnSet} from text file
#' 
#' @param data The name of a tab-separated file that contains the data.
#' 
#' @param sample A dataframe describing the samples (in lines).
#' 
#' @param indExpData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{fData()} table of the \code{MSnSet} object.
#' 
#' @param keyId The indice of the column containing the ID of entities 
#' (peptides or proteins)
#' 
#' @param namesOrigin xxxxxxxxxxx
#' 
#' @param logTransform A boolean value to indicate if the data have to be
#' log-transformed (Default is FALSE)
#' 
#' @param forceNA A boolean value to indicate if the 0 and NaN values of
#' intensity have to be replaced by NA (Default is FALSE)
#' 
#' @param typeOfData A string that indicates whether the dataset is about
#' 
#' @param parentProtId For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' proteins and is used to generate adjacency matrix and process to aggregation.
#' 
#' @param processes xxxx
#' 
#' @param pipelineType xxxx
#' 
#' @param analysis xxx
#' 
#' @return An instance of class \code{Features}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' data.file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata2")
#' data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata2")
#' sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' indExpData <- c(56:61)
#' namesOrigin <- colnames(data)[43:48]
#' parentId <- 'Protein_group_IDs'
#' keyid <- 'Sequence'
#' ft <- createFeatures(data, sample,indExpData,  keyId = keyid, namesOrigin = namesOrigin, typeOfData = "peptide", parentProtId = parentId, forceNA=TRUE)
#' 
#' @import Features
#' @importFrom utils installed.packages
#' @import SummarizedExperiment
#' 
#' @export
#' 
createFeatures <- function(data,
                           sample,
                           indExpData,
                           keyId=NULL,
                           namesOrigin = NULL,
                           logTransform=FALSE, 
                           forceNA=FALSE,
                           typeOfData,
                           parentProtId = NULL,
                           analysis,
                           processes = NULL,
                           pipelineType = NULL){
  
  if(missing(data))
    stop("'data' is missing.")
  if(missing(sample))
    stop("'sample' is missing.")
  if(missing(indExpData))
    stop("'indExpData' is missing.")
  if(missing(typeOfData))
    stop("'typeOfData' is missing.")
  
  
  
  
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
  SummarizedExperiment::colData(obj)@listData <- sample
  
  
  ## Replace all '.' by '_' in names
  #colnames(fd) <- gsub(".", "_", colnames(data)[indFData], fixed=TRUE)
  #colnames(Intensity) <- gsub(".", "_", colnames(data)[indExpData], fixed=TRUE)
  
  # As the function addOribingOfValues is based on the presence of NA in quanti data,
  # if forceNA is not set as TRUE, the previous function cannot ben run
  origin <- MultiAssayExperiment::DataFrame()
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
  
  
  daparVersion <- if (is.na(utils::installed.packages()["DAPAR2"])) 'NA' else utils::installed.packages()["DAPAR2",'Version']
  ProstarVersion <-if (is.na(utils::installed.packages()["Prostar2"])) 'NA' else utils::installed.packages()["Prostar2",'Version']
  
 
  metadata(obj) <- list(versions = list(Prostar_Version = ProstarVersion,
                                        DAPAR_Version = daparVersion),
                        parentProtId = parentProtId,
                        keyId = keyId,
                        params = list(),
                        typeOfData = typeOfData,
                        RawPValues = FALSE,
                        OriginOfValues = colnames(origin),
                        analysis = analysis,
                        pipelineType = pipelineType,
                        processes=c('original',processes)
  )
 
  
  return(obj)
}

