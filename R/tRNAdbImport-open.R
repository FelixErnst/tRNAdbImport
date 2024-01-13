#' @include tRNAdbImport.R
NULL

#' @name open_tdbID
#' @aliases open_tdbID open_mtdbID
#' 
#' @title Open a tRNA db entry in a browser
#' 
#' @description 
#' \code{open.tdbID} is a wrapper for \code{browseURL} and opens a tab for a
#' tRNAdb entry in a browser. Please note, that the tRNAdb server does not show
#' the entry right away without a session ID. open twice upon first use.
#'
#' @param tdbID a tRNA db
#' @param mtdbID a mtRNA db
#' @param dbURL the URL for the tRNAdb
#' 
#' @return opens a window in a default browser for tRNAdb entry selected
#' 
#' @importFrom utils browseURL
#' @importFrom httr2 url_build url_parse
#'
#' @export
#'
#' @examples
#' if(interactive()){
#'   open_tdbID("tdbD00000785")
#'   open_mtdbID("mtdbD00000907")
#' }
open_tdbID <- function(tdbID, dbURL = TRNA_DB_URL){
  utils::browseURL(
    httr2::url_build(
      httr2::url_parse(
        paste0(dbURL, "DataOutput/Result?ID=",
               tdbID))))
}
#' @rdname open_tdbID
#' @export
open_mtdbID <- function(mtdbID, dbURL = TRNA_DB_URL_MT){
  utils::browseURL(
    httr2::url_build(
      httr2::url_parse(
        paste0(dbURL, "DataOutput/Result?ID=",
               mtdbID))))
}
