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
#' @param gr the \code{GRanges} object to test
#' 
#' @return a logical value
#' 
#' @examples 
#' \donttest{
#' istRNAdbGRanges(gr)
#' }
NULL
#' @rdname istRNAdbGRanges
#' @export
setMethod(
  f = "istRNAdbGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr) .check_trnadb_granges(gr,
                                                    TRNADB_FEATURES))

# checks whether a GRanges object is tRNAdb compatible
.check_trnadb_granges <- function(gr,features){
  if(class(gr) != "GRanges"){
    stop("Input is not a GRanges object.",
         call. = FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. Please refer to the vignette of tRNAdbImport for ",
         "an exmaple on what information is expected.",
         call. = FALSE)
  }
  return(TRUE)
}