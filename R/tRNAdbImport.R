#' @title 
#' tRNAdbImport: Importing from to tRNAdb and mitotRNAdb as GRanges
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' title
#'
#' @references 
#' 
#' @docType package
#' @name tRNAdbImport
NULL


#' @import methods
#' @import assertive
#' @import httr
#' @import xml2
#' @import GenomicRanges
#' @import tRNAscanImport
NULL
requireNamespace("assertive")


# constants --------------------------------------------------------------------

#' @rdname tRNAdbImport
#' @export
TRNA_DB_URL <- "http://trna.bioinf.uni-leipzig.de"
#' @rdname tRNAdbImport
#' @export
TRNA_DB_URL_MT <- "http://mttrna.bioinf.uni-leipzig.de/mt"

TRNA_DB_VERIFIED <- c("verified sequence" = TRUE,
                      "unverified sequence" = FALSE,
                      "not verifiable sequence" = NA)
TRNA_DB_TYPE <- c("RNA","DNA")
TRNA_DB_ORIGIN <- c("plastid" = "chloro",
                    "mitochondrial" = "mito",
                    "allothers" = "allothers")