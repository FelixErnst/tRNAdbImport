#' @include tRNAdbImport.R
NULL

#' @name istRNAdbGRanges
#' @aliases istRNAdbGRanges
#' 
#' @title tRNAdb compatibility check
#' 
#' @description 
#' \code{istRNAdbGRanges} checks whether a GRanges object contains the 
#' information expected for a tRNAdb result.
#' 
#' @param x the \code{GRanges} object to test
#' 
#' @return a logical value
#' 
#' @examples 
#' gr <- import.tRNAdb(organism = "Saccharomyces cerevisiae",
#'               aminoacids = c("Phe","Ala"),
#'               anticodons = c("GAA"))
#' istRNAdbGRanges(gr)
NULL
#' @rdname istRNAdbGRanges
#' @export
setMethod(
  f = "istRNAdbGRanges",
  signature = signature(x = "GRanges"),
  definition = function(x) .check_trnadb_granges(x,TRNADB_FEATURES))

# checks whether a GRanges object is tRNAdb compatible
.check_trnadb_granges <- function(gr,features){
  if(!is(gr,"GRanges")){
    stop("Input is not a GRanges object.",
         call. = FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. The following columns are expected:\n'",
         paste(features, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(TRUE)
}
