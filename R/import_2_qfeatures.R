#' @title Creates an object of class \code{QFeatures} from text file.
#' 
#' @description 
#' 
#' Creates an object of class \code{QFeatures} from a
#' single tabulated-like file for quantitative and meta-data and a dataframe
#' for the samples description.
#'
#' @param data The name of a tab-separated file that contains the data.
#'
#' @param file `character(1)`. The name of a file xxx
#'
#' @param sample A dataframe describing the samples (in lines).
#'
#' @param indQData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{rowData()} table of the \code{QFeatures} object.
#'
#' @param keyId A `character(1)` or `numeric(1)` which is the indice of the column containing the ID of entities
#' (peptides or proteins)
#'
#' @param indQMetadata xxxxxxxxxxx
#'
#' @param force.na A `boolean` that indicates if the '0' and 'NaN' values of
#' quantitative values  must be replaced by 'NA' (Default is FALSE)
#'
#' @param typeDataset A string that indicates whether the dataset is about
#'
#' @param parentProtId A `character(1)` For peptide entities, a string which is the name of a column in rowData. It contains the id of parent
#' proteins and is used to generate adjacency matrix and process to aggregation.
#'
#' @param processes A vector of A `character()` which contains the name of processes
#' which has already been run on the data. Default is 'original'.
#'
#' @param typePipeline A `character(1)` The type of pipeline used with this dataset. The list of predefined
#' pipelines in DaparToolshed can be obtained with the function \code{pipelines()}. Default value is NULL
#'
#' @param analysis A `character(1)` which is the name of the MS study. 
#'
#' @param software A `character(1)` 
#' 
#' @param name A `character(1)` which is the name of the assay in the QFeatures object. 
#' Default is 'original'
#'
#' @return An instance of class \code{QFeatures}.
#'
#' @author Samuel Wieczorek
#'
#' @examples {
#' data.file <- system.file("extdata", "Exp1_R25_pept.txt", package="DaparToolshedData")
#' data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package="DaparToolshedData")
#' sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' ft <- createQFeatures(data = data, sample = sample, indQData = 56:61, keyId = 'Sequence', analysis = 'test',
#' indQMetadata = 43:48, typeDataset = 'peptide',
#' parentProtId = 'Protein_group_IDs', forceNA = TRUE, software = 'maxquant')
#' }
#'
#' @import QFeatures
#' @importFrom utils installed.packages
#' @import SummarizedExperiment
#'
#' @export
#' 
#' @rdname import-export-QFeatures
#'
createQFeatures <- function(data = NULL,
                            file = NULL,
                            sample,
                            indQData,
                            keyId = 'AutoID',
                            indQMetadata = NULL,
                            force.na = TRUE,
                            typeDataset,
                            parentProtId = NULL,
                            analysis = 'foo',
                            processes = NULL,
                            typePipeline = NULL,
                            software = NULL,
                            name = 'original'){


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
  
  if (is.numeric(indQData))
    indQData <- colnames(data)[indQData]
  
  if (is.numeric(indQMetadata))
    indQMetadata <- colnames(data)[indQMetadata]

  
  keyId <- ReplaceSpecialChars(keyId)
  typeDataset <- ReplaceSpecialChars(typeDataset)
  parentProtId <- ReplaceSpecialChars(parentProtId)
  analysis <- ReplaceSpecialChars(analysis)
  processes <- ReplaceSpecialChars(processes)
  typePipeline <- ReplaceSpecialChars(typePipeline)
  software <- ReplaceSpecialChars(software)


  
  
  if (keyId == 'AutoID'){
    auto <- rep(paste(typeDataset, "_", seq_len(nrow(data)), sep=""))
    data <- cbind(data, AutoID = auto)
    rownames(data) <- auto
  } else {
    rownames(data) <- data[,keyId]
  }
  # Creates the QFeatures object
  obj <- QFeatures::readQFeatures(data,
                                  ecol = indQData,
                                  name = 'original',
                                  fnames = keyId)
  
  ## Encoding the sample data
  sample <- lapply(sample, function(x){ ReplaceSpecialChars(x)})
  SummarizedExperiment::colData(obj)@listData <- sample
  
  
  # Get the quantitative metadata 
  tmp.qMetadata <- NULL
  if (!is.null(indQMetadata)){
    tmp.qMetadata <- data[, indQMetadata]
    tmp.qMetadata <- apply(tmp.qMetadata, 2, tolower)
    tmp.qMetadata <- apply(tmp.qMetadata, 
                           2,
                           function(x) 
                             gsub("\\s", "", x)
                           )
    tmp.qMetadata <- as.data.frame(tmp.qMetadata,
                                   stringsAsFactors = FALSE)
  }

  #browser()
  qMetadata <- BuildqMetadata(from = software,
                              level = typeDataset,
                              qdata = assay(obj),
                              conds = colData(obj)$Condition,
                              df = tmp.qMetadata)

  
  
  # Remove the identification columns which became useless
  rowData(obj[[1]])<- rowData(obj[[1]])[, -match(indQMetadata, colnames(rowData(obj[[1]])))]
  
  
  # Add the quantitative cell metadata info
  qMetadata(obj[['original']]) <- qMetadata
  

  if (force.na)
    obj <- zeroIsNA(obj, seq_along(obj))


  # Enrich the metadata for whole QFeatures object
  metadata(obj)$versions <- ProstarVersions()
  metadata(obj)$analysis <- list(analysis = analysis,
                                 typePipeline = typePipeline,
                                 processes = c('original', processes)
  )
  
  # Fill the metadata for the first assay
  typeDataset(obj[['original']]) <- typeDataset
  parentProtId(obj[['original']]) <- parentProtId
  idcol(obj[['original']]) <- keyId

  if (tolower(typeDataset) == 'peptide'){
    X<- makeAdjacencyMatrix(rowData(obj[[1]])[ ,parentProtId])
    rownames(X) <- rownames(rowData(obj[[1]]))
    adjacencyMatrix(obj[[1]]) <- X
  }

return(obj)
}

