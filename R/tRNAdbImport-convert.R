#' @include tRNAdbImport.R
NULL

#' @name tRNAdb2GFF
#' 
#' 
#' @title Convert result to GFF format compatible result
#' 
#' @description
#' To facilitate interaction with standard \code{GRanges} objects and repective
#' functions and ecosystem available, \code{tRNAdb2GFF} convert to standard
#' \code{GFF} format returning a \code{GRanges} object.
#' 
#' @param input a GRanges object which passes the \code{istRNAdbGRanges} check
#'
#' @return a GRanges object containing the information from the tRNA db converted
#' and santized for GFF compatbility.
#' 
#' @export
#' 
#' @examples
#' tRNAs <- import.tRNAdb(organism = "Saccharomyces cerevisiae",
#'                        aminoacids = c("Phe","Ala"))
#' tRNAdb2GFF(tRNAs)
tRNAdb2GFF <- function(input) {
  .check_trnadb_granges(input, TRNADB_FEATURES)
  tRNAdb <- input
  # patch GRanges object with necessary columns for gff3 comptability
  S4Vectors::mcols(tRNAdb)$tRNA_seq <- 
    as.character(S4Vectors::mcols(tRNAdb)$tRNA_seq)
  S4Vectors::mcols(tRNAdb)$tRNA_str <- 
    as.character(S4Vectors::mcols(tRNAdb)$tRNA_str)
  S4Vectors::mcols(tRNAdb)$ID <- S4Vectors::mcols(tRNAdb)$tRNAdb_ID
  S4Vectors::mcols(tRNAdb)$type <- "tRNA"
  S4Vectors::mcols(tRNAdb)$type <- 
    as.factor(S4Vectors::mcols(tRNAdb)$type)
  S4Vectors::mcols(tRNAdb)$source <- "tRNAdb"
  S4Vectors::mcols(tRNAdb)$source <- 
    as.factor(S4Vectors::mcols(tRNAdb)$source)
  S4Vectors::mcols(tRNAdb)$score <- NA
  S4Vectors::mcols(tRNAdb)$score <- 
    as.numeric(S4Vectors::mcols(tRNAdb)$score)
  S4Vectors::mcols(tRNAdb)$phase <- NA
  S4Vectors::mcols(tRNAdb)$phase <- 
    as.integer(S4Vectors::mcols(tRNAdb)$phase)
  S4Vectors::mcols(tRNAdb)$score <- 
    as.integer(S4Vectors::mcols(tRNAdb)$phase)
  # arrange columns in correct order
  S4Vectors::mcols(tRNAdb) <- 
    cbind(S4Vectors::mcols(tRNAdb)[,c("source",
                                      "type",
                                      "score",
                                      "phase",
                                      "ID")],
          S4Vectors::mcols(tRNAdb)[,-which(colnames(
            S4Vectors::mcols(tRNAdb)) %in% 
              c("source",
                "type",
                "score",
                "phase",
                "ID"))])
  return(tRNAdb)
}
