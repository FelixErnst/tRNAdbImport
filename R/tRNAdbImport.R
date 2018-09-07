#' @title 
#' tRNAdbImport: Importing from to tRNAdb and mitotRNAdb as GRanges
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' title
#'
#' @references 
#' Jühling F, Mörl M, Hartmann RK, Sprinzl M, Stadler PF, Pütz J. 2009. "tRNAdb 
#' 2009: compilation of tRNA sequences and tRNA genes." Nucleic Acids Research, 
#' Volume 37 (suppl_1): D159–162. doi:10.1093/nar/gkn772.
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