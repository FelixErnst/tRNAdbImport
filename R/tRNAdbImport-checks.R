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
  definition = function(x) .check_trnadb_granges(x, TRNADB_FEATURES))

# checks whether a GRanges object is tRNAdb compatible
.check_trnadb_granges <- function(gr,features){
  if(!is(gr,"GRanges")){
    warning("Input is not a GRanges object.", call. = FALSE)
    return(FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    warning("Input GRanges object does not meet the requirements of the ",
            "function. The following columns are expected:\n'",
            paste(features, collapse = "', '"),
            "'.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#input type checks
.checkValueValidity <- function(value, checkValues,
                                .xvalue = .get_name_in_parent(value)){
  if(!all(value %in% checkValues)){
    stop("'",gsub("\"","",.xvalue),
         "' must be one of the following values: '",
         paste(checkValues, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}

.get_name_in_parent <- function(x) {
  .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
  paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}
