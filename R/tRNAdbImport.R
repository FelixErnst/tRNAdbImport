#' @title 
#' tRNAdbImport: Importing from to tRNAdb and mitotRNAdb as GRanges
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' The tRNAdb and mttRNAdb (Jühling et al. 2009) is a compilation of tRNA 
#' sequences and tRNA genes. It is a follow up version of the database of 
#' Sprinzl et al. 2005.
#' 
#' Using `tRNAdbImport` the tRNAdb can be accessed as outlined on the website 
#' [http://trna.bioinf.uni-leipzig.de/](http://trna.bioinf.uni-leipzig.de/) and 
#' the results are returned as a `GRanges` object.
#'
#' @section Manual:
#' Please refer to the tRNAdbImport vignette for an example how to work and 
#' use the package: \href{../doc/tRNAdbImport.html}{tRNAdbImport} 
#' 
#' @seealso [import.tRNAdb()] for examples
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
#' @import Biostrings
#' @import Structstrings
#' @import Modstrings
#' @import tRNA
NULL

# constants --------------------------------------------------------------------

TRNA_DB_VERIFIED <- c("verified sequence" = TRUE,
                      "unverified sequence" = FALSE,
                      "not verifiable sequence" = NA)
TRNA_DB_TYPE <- c("RNA","DNA")
TRNA_DB_ORIGIN <- c("plastid" = "chloro",
                    "mitochondrial" = "mito",
                    "allothers" = "allothers")

TRNADB_FEATURES <- c(
  tRNA:::TRNA_FEATURES,
  "tRNAdb_ID",
  "tRNAdb",
  "tRNAdb_organism",
  "tRNAdb_strain",
  "tRNAdb_taxonomyID",
  "tRNAdb_verified"
)
