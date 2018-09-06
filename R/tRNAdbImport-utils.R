#' @include tRNAdbImport.R
NULL

# assemble URL for accessing tRNA db
.trnadb_url <- function(url, path = "") {
  url <- httr::modify_url(url, path = path)
  return(url)
}

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