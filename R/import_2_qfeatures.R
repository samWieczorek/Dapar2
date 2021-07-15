
#' @title Standardize names
#' 
#' @description Replace some characters in names by 'underscore'
#' 
#' @param x A vector of strings to be processed
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
ReplaceSpecialChars <- function(x){
  if (is.null(x))
    return(x)

  val <- x
  for (c in c(".", ' ', '-'))
    val <- gsub(c, '_', x, fixed=TRUE)
  val
}
  

#' Builds an object of class \code{QFeatures} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description.
#' 
#' @title Creates an object of class \code{QFeatures} from text file
#' 
#' @param data The name of a tab-separated file that contains the data.
#' 
#' @param file The name of a file xxx
#' 
#' @param sample A dataframe describing the samples (in lines).
#' 
#' @param indQData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{rowData()} table of the \code{QFeatures} object.
#' 
#' @param keyId The indice of the column containing the ID of entities 
#' (peptides or proteins)
#' 
#' @param indQMetadata xxxxxxxxxxx
#' 
#' @param logTransform A boolean value to indicate if the data have to be
#' log-transformed (Default is FALSE)
#' 
#' @param forceNA A boolean value to indicate if the 0 and NaN values of
#' intensity have to be replaced by NA (Default is FALSE)
#' 
#' @param typeDataset A string that indicates whether the dataset is about
#' 
#' @param parentProtId For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' proteins and is used to generate adjacency matrix and process to aggregation.
#' 
#' @param processes xxxx
#' 
#' @param typePipeline The type of pipeline used with this dataset. The list of predefined
#' pipelines in DaparToolshed can be obtained with the function \code{pipelines()}. Default value is NULL
#' 
#' @param analysis The name of the MS analysis
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
#' ft <- createQFeatures(data, sample, 
#' indQData = 56:61, 
#' keyId = 'Sequence', 
#' analysis = 'test,
#' logTransform = TRUE,
#' indQMetadata = 43:48, 
#' typeDataset = 'peptide', 
#' parentProtId = 'Protein_group_IDs', 
#' forceNA = TRUE, 
#' software = 'maxquant')
#' 
#' @import QFeatures
#' @importFrom utils installed.packages
#' @import SummarizedExperiment
#' 
#' @export
#' 
createQFeatures <- function(data,
                            file,
                            sample,
                            indQData,
                            keyId = 'AutoID',
                            indQMetadata = NULL,
                            logTransform = FALSE, 
                            forceNA = TRUE,
                            typeDataset,
                            parentProtId = NULL,
                            analysis = 'foo',
                            processes = NULL,
                            typePipeline = NULL,
                            software = NULL){
  
  
  #Check parameters validity
  if(missing(data) && missing(file))
    stop("Either 'data' or 'file' is required")
  else if (!missing(data) && !missing(file))
    stop("Only 'data' or 'file' is required at a time. Please choose one of them.")
  
  if (!missing(data) && class(data) != "data.frame")
    stop("'data' must be a data.frame")
  
  if (!missing(file)){
    if(class(file) != "xxx")
      stop("'file' must be a connection")
    else
      data <- read.table(file, 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE)
  }
  
  
  
  if(missing(sample))
    stop("'sample' is required")
  else if (class(sample) != "data.frame")
    stop("'sample' must be a data.frame")
  
  
  if(missing(indQData))
    stop("'indQData' is required")
  else if (!is.numeric(indQData))
    stop("'indQData' must be a vector of integer")
  
  if(missing(indQMetadata))
    stop("'indQMetadata' is required")
  else if (!is.numeric(indQMetadata))
    stop("'indQMetadata' must be a vector of integer")
  
  if (!is.null(keyId) && !is.character(keyId))
    stop("'keyId' must be either NULL nor a string")
  
  if(missing(typeDataset))
    stop("'typeDataset' is required")
  
  #Standardize all colnames
  colnames(data) <- ReplaceSpecialChars(colnames(data))
  
  keyId <- ReplaceSpecialChars(keyId)
  typeDataset <- ReplaceSpecialChars(typeDataset)
  parentProtId <- ReplaceSpecialChars(parentProtId)
  analysis <- ReplaceSpecialChars(analysis)
  processes <- ReplaceSpecialChars(processes)
  typePipeline <- ReplaceSpecialChars(typePipeline)
  software <- ReplaceSpecialChars(software)
  
  
  if (keyId == 'AutoID')
    data <- cbind(data, 
                  AutoID = rep(paste(typeDataset, "_", 1:nrow(data), sep=""))
    )
  
  obj <- QFeatures::readQFeatures(data, 
                                  ecol = indQData, 
                                  name = 'original',
                                  fnames = keyId)
  
  ## Encoding the sample data
  sample <- lapply(sample, function(x){ ReplaceSpecialChars(x)})
  SummarizedExperiment::colData(obj)@listData <- sample
  
  
  # Get the metacell info
  tmp.metacell <- NULL
  if (!is.null(indQMetadata)){
    tmp.metacell <- data[, indQMetadata]
    tmp.metacell <- apply(tmp.metacell, 2, tolower)
    tmp.metacell <- as.data.frame(apply(tmp.metacell, 2, 
                                        function(x) gsub("\\s", "", x)),
                                  stringsAsFactors = FALSE)
  }
  
  #browser()
  metacell <- BuildMetaCell(from = software,
                            level = typeDataset,
                            qdata = assay(obj), 
                            conds = colData(obj)$Condition,
                            df = tmp.metacell)
  
  
  # Add the quantitative cell metadata info
  rowData(obj[['original']]) <- cbind(rowData(obj[['original']]), 
                                      as.data.frame(metacell, rownames = colnames(data)[indQMetadata])
  )
  
  
  if (isTRUE(forceNA))
    obj <- zeroIsNA(obj, seq_along(obj))
  
  
  # Fill the metadata for whole object
  metadata(obj)$versions <- GetProstarVersions()
  metadata(obj)$keyId <- keyId
  metadata(obj)$params <- list()
  metadata(obj)$RawPValues <- FALSE
  metadata(obj)$qMetadata_names <- colnames(data)[indQMetadata]
  metadata(obj)$analysis <- analysis
  metadata(obj)$typePipeline <- typePipeline
  metadata(obj)$processes <- c('original', processes)
  
  # Fill the metadata for the first assay
  metadata(obj[['original']])$typeDataset <- typeDataset
  metadata(obj[['original']])$parentProtId <- parentProtId
  
  
  if (tolower(typeDataset) == 'peptide'){
    obj <- SetAdjMat(obj, 1)
    obj <- SetConnectedComps(obj, 1)
  }
  
  
  if (isTRUE(logTransform)) {
    obj <- addAssay(obj, 
                    logTransform(obj[['original']]),
                    name = "original_log")
    obj <- addAssayLinkOneToOne(obj, 
                                from = "original",
                                to = "original_log")
  }
  
  
  return(obj)
}



#' 
#' 
#' #' @title Creates an object of class \code{QFeatures} from an object of class \code{MSnSet}
#' #' 
#' #' @description xxxx
#' #' 
#' #' @param obj xxx.
#' #' 
#' #' @param analysis xxx
#' #' 
#' #' @param parentProtId For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' #' proteins and is used to generate adjacency matrix and process to aggregation.
#' #' 
#' #' @param keyId The indice of the column containing the ID of entities 
#' #' (peptides or proteins)
#' #' 
#' #' @param pipelineType xxxx
#' #' 
#' #' @param processes xxxx
#' #' 
#' #' @return An instance of class \code{QFeatures}.
#' #' 
#' #' @author Samuel Wieczorek
#' #' 
#' #' @examples 
#' #' library(QFeatures)
#' #' library(MSnbase)
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' obj <- Exp1_R25_pept
#' #' parentId <- 'Protein_group_IDs'
#' #' keyid <- 'Sequence'
#' #' ft <- convertMSnset2QFeatures(obj, 'conv',parentId, keyid )
#' #' 
#' #' @importFrom Biobase exprs fData pData
#' #' 
#' #' @importFrom QFeatures readQFeatures
#' #' 
#' #' @export
#' #' 
#' convertMSnset2QFeatures <- function(obj, analysis, parentProtId, keyId, pipelineType = NULL, processes = NULL) {
#'   
#'   
#'   if (class(obj) != 'MSnSet')
#'     stop("This dataset is not a MSnset file.")
#'   
#'   if(missing(analysis))
#'     stop("'analysis' is required.")
#'   if(missing(parentProtId))
#'     stop("'parentProtId' is required.")
#'   if(missing(keyId))
#'     stop("'keyId' is required.")
#'   # if(missing(pipelineType))
#'   #   stop("'pipelineType' is required.")
#'   # if(missing(processes))
#'   #   stop("'processes' is required.")
#'   
#'   
#'   df <- cbind(Biobase::fData(obj), Biobase::exprs(obj))
#'   i <- (ncol(Biobase::fData(obj))+1):(ncol(df))
#'   feat <- QFeatures::readQFeatures(df, ecol = i, sep = "\t", name = "original", fnames = keyId)
#'   feat[['original']] <- SetTypeDataset(feat[['original']], obj@experimentData@other$typeOfData)
#'   
#'   ## Encoding the sample data
#'   sample <- lapply(Biobase::pData(obj),function(x){ gsub(".", "_", x, fixed=TRUE)})
#'   SummarizedExperiment::colData(feat)@listData <- sample
#'   
#'   feat <- QFeatures::zeroIsNA(feat,seq_along(feat))
#'   
#'   if (is.null(obj@experimentData@other$OriginOfValues))
#'   {
#'     warning("The MSnset file odes not contain any information about the origin of values.")
#'     warning("So, this convert tool cannot be used. Please use the convert GUI from raw datasets in Prostar.")
#'     return(NULL)
#'   } else {
#'     origin <- obj@experimentData@other$OriginOfValues
#'   }
#'   
#'   daparVersion <- if (is.na(utils::installed.packages()["DAPAR2"])) 'NA' else utils::installed.packages()["DAPAR2",'Version']
#'   ProstarVersion <-if (is.na(utils::installed.packages()["Prostar2"])) 'NA' else utils::installed.packages()["Prostar2",'Version']
#'   
#'   
#'   metadata(feat) <- list(versions = list(Prostar_Version = ProstarVersion,
#'                                          DAPAR_Version = daparVersion),
#'                          parentProtId = parentProtId,
#'                          keyId = keyId,
#'                          params = list(),
#'                          RawPValues = obj@experimentData@other$RawPValues,
#'                          OriginOfValues = origin,
#'                          analysis = analysis,
#'                          pipelineType = pipelineType,
#'                          processes=c('original',processes)
#'   )
#'   
#'   return(feat)
#'   
#'   
#' }
