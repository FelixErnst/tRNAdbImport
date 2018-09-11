#' @include tRNAdbImport.R
NULL

#' @name open_tdbID
#' @aliases open_tdbID open_mtdbID
#' 
#' @title Open a tRNA db entry in a browser
#' 
#' @description 
#' \code{open.tdbID} is a wrapper for \code{browseURL} and opens a tab for a
#' tRNAdb entry in a browser.
#'
#' @param tdbID a tRNA db
#' @param mtdbID a mtRNA db
#' @param dbURL the URL for the tRNAdb
#'
#' @export
#'
#' @examples
#' \donttest{
#' open_tdbID("tdbD00000785")
#' open_mtdbID("mtdbD00000907")
#' }
open_tdbID <- function(tdbID,
                       dbURL = TRNA_DB_URL){
  browseURL(paste0(httr::modify_url(dbURL),
                   "DataOutput/Result?ID=",
                   tdbID))
}
#' @rdname open_tdbID
#' @export
open_mtdbID <- function(mtdbID,
                        dbURL = TRNA_DB_URL_MT){
  browseURL(paste0(httr::modify_url(dbURL),
                   "DataOutput/Result?ID=",
                   mtdbID))
}