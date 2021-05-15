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
setMEC <- function(origin, qData, conds=NULL){
  
  if (missing(origin))
    stop("'origin' is required.")
  
  if (missing(qData))
    stop("'qData' is required.")
  if (missing(conds))
    stop("'conds' is required.")
  
  
  if (is.null(conds)){
    warning("'conds' is not set. The object is unchanged.")
    return(origin)
  }
  
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
#' @param obj An object of class \code{QFeatures}
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
#' origin <- c("Identification_type_C_R1", "Identification_type_C_R2", 
#' "Identification_type_C_R3", "Identification_type_D_R1", 
#' "Identification_type_D_R2", "Identification_type_D_R3")
#' obj <- addOriginOfValues(Exp1_R25_pept, 2, origin)
#' }
#' 
#' @importFrom methods as
#' 
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' 
#' @export
#' 
addOriginOfValues <- function(obj, i, namesOrigin = NULL){
  
  if (!is.null(namesOrigin))
  {
    if ( ncol(assay(obj,i)) != length(namesOrigin)){
      stop("The number of samples in the assay must be equal to the length of 'namesOrigin'")
    }
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
                           qData = assay(obj[[i]]), 
                           conds = SummarizedExperiment::colData(obj)$Condition)
  
  
  return(OriginOfValues)
}






#' Builds an object of class \code{QFeatures} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description.
#' 
#' @title Creates an object of class \code{QFeatures} from text file
#' 
#' @param data The name of a tab-separated file that contains the data.
#' 
#' @param sample A dataframe describing the samples (in lines).
#' 
#' @param indExpData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{rowData()} table of the \code{QFeatures} object.
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
#' @return An instance of class \code{QFeatures}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' data.file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata2")
#' data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata2")
#' sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' indExpData <- 56:61
#' parentId <- 'Protein_group_IDs'
#' keyId <- 'Sequence'
#' analysis <- 'test'
#' indexForMetacell <- c(43:48)
#' ft <- createQFeatures(data, sample, indExpData, keyId = keyid, analysis=analysis, logTransform = TRUE,
#' indexForMetacell = indexForMetacell, typeOfData = "peptide", parentProtId = parentId, forceNA=TRUE, software = 'maxquant')
#' 
#' @import QFeatures
#' @importFrom utils installed.packages
#' @import SummarizedExperiment
#' 
#' @export
#' 
createQFeatures <- function(data,
                            sample,
                            indExpData,
                            keyId=NULL,
                            indexForMetacell = NULL,
                            logTransform=FALSE, 
                            forceNA=FALSE,
                            typeOfData,
                            parentProtId = NULL,
                            analysis='foo',
                            processes = NULL,
                            pipelineType = NULL,
                            software = NULL){
  
  if(missing(data))
    stop("'data' is required")
  if(missing(sample))
    stop("'sample' is required")
  if(missing(indExpData))
    stop("'indExpData' is required")
  if(missing(typeOfData))
    stop("'typeOfData' is required")
  
  
  colnames(data) <- gsub(".", "_", colnames(data), fixed=TRUE)
  keyId <- gsub(".", "_", keyId, fixed=TRUE)
  typeOfData <- gsub(".", "_", typeOfData, fixed=TRUE)
  parentProtId <- gsub(".", "_", parentProtId, fixed=TRUE)
  analysis <- gsub(".", "_", analysis, fixed=TRUE)
  processes <- gsub(".", "_", processes, fixed=TRUE)
  pipelineType <- gsub(".", "_", pipelineType, fixed=TRUE)
  software <- gsub(".", "_", software, fixed=TRUE)
  
  if (is.null(keyId) || keyId == '' || nchar(keyId)==0 || keyId == 'AutoID') {
    obj <- QFeatures::readQFeatures(data, 
                                    ecol=indExpData, 
                                    name='original')
  } else {
    obj <- QFeatures::readQFeatures(data, 
                                    ecol=indExpData, 
                                    name='original',
                                    fnames = keyId)
  }
  
  
  ## Encoding the sample data
  sample <- lapply(sample, function(x){ gsub(".", "_", x, fixed=TRUE)})
  SummarizedExperiment::colData(obj)@listData <- sample
  
  
  # Get the metacell info
  tmp.metacell <- NULL
  if (!is.null(indexForMetacell)){
    tmp.metacell <- data[, indexForMetacell]
    tmp.metacell <- apply(tmp.metacell,2,tolower)
    tmp.metacell <- as.data.frame(apply(tmp.metacell, 2, function(x) gsub("\\s", "", x)),
                              stringsAsFactors = FALSE)
  }
  
  #browser()
  metacell <- BuildMetaCell(from = software,
                            level = typeOfData,
                            qdata = assay(obj), 
                            conds = colData(obj)$Condition,
                            df = tmp.metacell)

  

  rowData(obj[['original']]) <- cbind(rowData(obj[['original']]), metacell)
  
  obj <- zeroIsNA(obj, seq_along(obj))
  
  
  daparVersion <- if (is.na(utils::installed.packages()["DAPAR2"])) 'NA' else utils::installed.packages()["DAPAR2",'Version']
  ProstarVersion <-if (is.na(utils::installed.packages()["Prostar2"])) 'NA' else utils::installed.packages()["Prostar2",'Version']
  
  metadata(obj) <- list(versions = list(Prostar_Version = ProstarVersion,
                                        DAPAR_Version = daparVersion),
                        parentProtId = parentProtId,
                        keyId = keyId,
                        params = list(),
                        RawPValues = FALSE,
                        names_metacell = colnames(metacell),
                        analysis = analysis,
                        pipelineType = pipelineType,
                        processes=c('original',processes)
  )
  
  metadata(obj[['original']])$typeOfData <- typeOfData
  
  
  if (tolower(typeOfData) == 'peptide')
  {
    obj <- addListAdjacencyMatrices(obj, 1)
    obj <- addConnexComp(obj, 1)
  }
  
  
  if (isTRUE(logTransform)) {
    obj <- addAssay(obj, logTransform(obj[['original']]),name = "original_log")
    obj <- addAssayLinkOneToOne(obj, from = "original", to = "original_log")
  }
  
  return(obj)
}





#' @title Creates an object of class \code{QFeatures} from an object of class \code{MSnSet}
#' 
#' @description xxxx
#' 
#' @param obj xxx.
#' 
#' @param analysis xxx
#' 
#' @param parentProtId For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' proteins and is used to generate adjacency matrix and process to aggregation.
#' 
#' @param keyId The indice of the column containing the ID of entities 
#' (peptides or proteins)
#' 
#' @param pipelineType xxxx
#' 
#' @param processes xxxx
#' 
#' @return An instance of class \code{QFeatures}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' library(QFeatures)
#' library(MSnbase)
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' parentId <- 'Protein_group_IDs'
#' keyid <- 'Sequence'
#' ft <- convertMSnset2QFeatures(obj, 'conv',parentId, keyid )
#' 
#' @importFrom Biobase exprs fData pData
#' 
#' @importFrom QFeatures readQFeatures
#' 
#' @export
#' 
convertMSnset2QFeatures <- function(obj, analysis, parentProtId, keyId, pipelineType = NULL, processes = NULL) {
  
  
  if (class(obj) != 'MSnSet')
    stop("This dataset is not a MSnset file.")
  
  if(missing(analysis))
    stop("'analysis' is required.")
  if(missing(parentProtId))
    stop("'parentProtId' is required.")
  if(missing(keyId))
    stop("'keyId' is required.")
  # if(missing(pipelineType))
  #   stop("'pipelineType' is required.")
  # if(missing(processes))
  #   stop("'processes' is required.")
  
  
  df <- cbind(Biobase::fData(obj), Biobase::exprs(obj))
  i <- (ncol(Biobase::fData(obj))+1):(ncol(df))
  feat <- QFeatures::readQFeatures(df, ecol = i, sep = "\t", name = "original", fnames = keyId)
  metadata(feat[['original']])$typeOfData <- obj@experimentData@other$typeOfData
  
  
  ## Encoding the sample data
  sample <- lapply(Biobase::pData(obj),function(x){ gsub(".", "_", x, fixed=TRUE)})
  SummarizedExperiment::colData(feat)@listData <- sample
  
  feat <- QFeatures::zeroIsNA(feat,seq_along(feat))
  
  if (is.null(obj@experimentData@other$OriginOfValues))
  {
    warning("The MSnset file odes not contain any information about the origin of values.")
    warning("So, this convert tool cannot be used. Please use the convert GUI from raw datasets in Prostar.")
    return(NULL)
  } else {
    origin <- obj@experimentData@other$OriginOfValues
  }
  
  daparVersion <- if (is.na(utils::installed.packages()["DAPAR2"])) 'NA' else utils::installed.packages()["DAPAR2",'Version']
  ProstarVersion <-if (is.na(utils::installed.packages()["Prostar2"])) 'NA' else utils::installed.packages()["Prostar2",'Version']
  
  
  metadata(feat) <- list(versions = list(Prostar_Version = ProstarVersion,
                                         DAPAR_Version = daparVersion),
                         parentProtId = parentProtId,
                         keyId = keyId,
                         params = list(),
                         RawPValues = obj@experimentData@other$RawPValues,
                         OriginOfValues = origin,
                         analysis = analysis,
                         pipelineType = pipelineType,
                         processes=c('original',processes)
  )
  
  return(feat)
  
  
}
