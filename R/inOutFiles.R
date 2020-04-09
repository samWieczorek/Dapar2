

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
#' exprsFile <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata = utils::read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE)
#' indExpData <- c(56:61)
#' indFData <- c(1:55,62:71)
#' keyid <- 'Sequence'
#' createMSnset(exprsFile, metadata,indExpData,  indFData, keyid, indexForOriginOfValue = NULL, typeOfData = "peptide")
#' @importFrom MSnbase MSnSet
#' @importFrom utils read.table
#' @export


