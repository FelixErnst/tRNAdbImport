#' @include tRNAdbImport.R
NULL

#input type checks
.checkValueValidity <- function(value,
                                checkValues,
                                .xvalue = assertive::get_name_in_parent(value)){
  if(!all(value %in% checkValues)){
    stop("'",.xvalue,
         "' must be one of the following values: '",
         paste(checkValues, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}