#' @include tRNAdbImport.R
NULL

#' @name import.tRNAdb.ce
#' @aliases import.tRNAdb.ce import.mttRNAdb.ce import.tRNAdb.ce.id import.mttRNAdb.ce.id
#' import.tRNAdb.ce.blast import.mttRNAdb.ce.blast tRNAdbce2GFF
#' 
#' @title Importing information from the tRNADB-CE as GRanges object
#' 
#' @description
#' For importing information from tRNAdb-CE multiple functions are available. 
#' They share a common syntax where possible. Specialized calls are available 
#' for id search and blast search.
#' 
#' The search type are fixed for \code{and} searches, only.
#'
#' @param organism a organism name as a character string 
#' @param strain a strain information as a character string
#' @param taxonomyID organism and strain information as a taxonom ID
#' @param aminoacids a character vector of amino acids as a three letter code 
#' @param anticodons a character vector of anticodon sequences 
#' @param sequences a named (1-15) list of sequences, which are used for the
#' search
#' @param structures a named (1-15) list of structures, which are used for the
#' search. Please use the \code{\(\)} or \code{><} dot bracket annotation.
#' @param reference a reference as a character string
#' @param comment a comment as a character string
#' @param pubmed a pubmed ID
#' @param genes a gene name as a character string
#' @param tdbID a tRNAdb ID 
#' @param mtdbID a mtRNAdb ID 
#' @param blastSeq a sequence to use for a blast search
#' @param database "RNA" or "DNA"
#' @param origin one ore more of "plastid", "mitochondrial" or "allothers"
#' @param dbURL the URL of the tRNA db
#' @param verbose whether to report verbose information from the httr2 calls
#'
#' @return a GRanges object containing the information from the tRNA db
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom stringr str_detect str_split str_locate_all
#' 
#' @export
#'
#' @examples
NULL

#' @rdname import.tRNAdb.ce
#' @export
TRNA_DB_CE_URL <- "https://trna.ie.niigata-u.ac.jp/"

# constants --------------------------------------------------------------------

TRNA_DB_CE_TYPE <- c("Bacteria", "Archaea", "WGS", "Plant", "Fungi", "Virus",
                     "Phage", "Plasmid", "Chloroplast", "ENV", "SRA")

# ------------------------------------------------------------------------------

.assemble_args_for_tRNA_db_ce_search <- function(keyword,
                                                 sequenceID,
                                                 genomeID, 
                                                 blastSequence,
                                                 eValue, 
                                                 maxResults,
                                                 blastOptions,
                                                 types,
                                                 reliable_targets){
  # base values
  args <- list()
  #
  if(!missing(sequenceID) && all(!is.na(sequenceID)) && all(sequenceID != "")){
    args[["KEY_WORD_SEQ"]] <- paste0(sequenceID, collapse = " ")
  }
  #
  if(!missing(genomeID) && all(!is.na(genomeID)) && all(genomeID != "")){
    args[["KEY_WORD_SEQ"]] <- paste0(genomeID, collapse = " ")
  }
  #
  if(!missing(reliable_targets) && length(reliable_targets) == 1L){
    args[["VTYPE"]] <- reliable_targets ? 1L : 0L
  }
  #
  if(!missing(keyword) && all(!is.na(keyword)) && all(keyword != "")){
    args[["KEY_WORD_NA"]] <- keyword
  }
  #
  if(!missing(blastSequence) && length(blastSequence) == 1L &&
     blastSequence != ""){
    args[["SEQ"]] <- blastSequence
  }
  #
  if(!missing(eValue) && length(eValue) == 1L){
    args[["EVAL"]] <- eValue
  }
  #
  if(!missing(maxResults) && length(maxResults) == 1L){
    args[["MAXVAL"]] <- maxResults
  }
  #
  if(!missing(blastOptions) && length(blastOptions) == 1L && blastOptions != ""){
    args[["BOP"]] <- blastOptions
  }
  #
  if(!missing(types) && all(!is.na(types)) && all(types != "")){
    args[["DYTPE"]] <- types
  }
  #
  args[["KEY_SEL_NA"]] <- 1L
  args[["SEARCH_JOIN"]] <- 1L
  args[["SEARCH_JOIN"]] <- 1L
  #
  args[["KEY_SEL_AA"]] <- 2L
  args[["KEY_SEL_AC"]] <- 2L
  args[["KEY_SEL_SEQ"]] <- 2L
  args[["KEY_SEL_GEN"]] <- 2L
  args[["ADV"]] <- 1L
  #
  args
}


#' @rdname import.tRNAdb.ce
#' @export
import.tRNAdb.ce.id <- function(sequenceID, genomeID,
                                types = c("Bacteria"), 
                                reliable_targets = TRUE,
                                dbURL = TRNA_DB_CE_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  if(!.is_non_empty_character(sequenceID)){
    stop("'sequenceID' must contain only non empty character values.")
  }
  if(!.is_non_empty_character(genomeID)){
    stop("'genomeID' must contain only non empty character values.")
  }
  if(!.is_a_bool(reliable_targets)){
    stop("'reliable_targets' must be TRUE or FALSE")
  }
  types <- match.arg(types,
                     c(TRNA_DB_CE_TYPE,NA_character_), 
                     several.ok = TRUE)
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_ce_search(sequenceID = sequenceID,
                                               genomeID = genomeID,
                                               types = types,
                                               reliable_targets = reliable_targets)
  args[["submit"]] <- "GO"
  args <- args[c("KEY_WORD_NA",
                 "KEY_SEL_NA",
                 "KEY_SEL_AA",
                 "KEY_SEL_AC",
                 "KEY_SEL_SEQ",
                 "KEY_SEL_GEN",
                 "ADV",
                 "VTYPE",
                 "DTYPE",
                 "submit")]
  # get result and return as GRanges
  res <- .get_trna_db_ce_result(args = args,
                                dbURL = dbURL,
                                dbFunction = "whole_seq_list_adv.cgi",
                                verbose = verbose)
  return(.convert_tRNAdb_ce_result_to_GRanges(res))
}

#' @rdname import.tRNAdb.ce
#' @export
import.tRNAdb.ce.blast <- function(blastSeq,
                                   eValue = 10,
                                   maxResults = 50,
                                   blastOptions = "",
                                   types = c("Bacteria"),
                                   dbURL = TRNA_DB_CE_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  if(.is_a_single_number(eValue)){
    stop("'eValue' must be single numeric value")
  }
  if(.is_a_single_number(maxResults)){
    stop("'maxResults' must be single numeric value")
  }
  blastSeq <- as.character(blastSeq) # in case it is a single DNAString*
  if(!.is_non_empty_string(blastSeq)){
    stop("'blastSeq' must be a single non empty character value.")
  }
  types <- match.arg(types,
                     c(TRNA_DB_CE_TYPE,NA_character_), 
                     several.ok = TRUE)
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_ce_search(blastSequence = blastSeq,
                                               eValue = eValue,
                                               maxResults = maxResults,
                                               blastOptions = blastOptions,
                                               types = types,
                                               reliable_targets = reliable_targets)
  args[["submit"]] <- "GO"
  args[["STYPE"]] <- "blastn"
  args <- args[c("SEQ",
                 "STYPE",
                 "EVAL",
                 "MAXVAL",
                 "BOP",
                 "DTYPE",
                 "submit")]
  # get result and return as GRanges
  res <- .get_trna_db_ce_result(args = args,
                                dbURL = dbURL,
                                dbFunction = "whole_seq_search.cgi",
                                verbose = verbose)
  return(.convert_tRNAdb_ce_result_to_GRanges(res))
}


.assemble_args_for_tRNA_db_ce_pattern <- function(upstream,
                                                  accstem,
                                                  np89,
                                                  dstem,
                                                  dloop,
                                                  dstem_rev,
                                                  np26,
                                                  acsteam,
                                                  acloop,
                                                  acstem_rev,
                                                  vloop,
                                                  tstem,
                                                  tloop,
                                                  tstem_rev,
                                                  accstem_rev,
                                                  np7376,
                                                  downstream,
                                                  accstem_size,
                                                  types){
  # base values
  args <- as.list(environment())
  args[["accstem_size"]] <- NULL
  args[["types"]] <- NULL
  #
  names(args) <- paste0("TSEQ_",seq(1,17))
  #
  if(!missing(types) && all(!is.na(types)) && all(types != "")){
    args[["DYTPE"]] <- types
  }
  #
  if(!missing(accstem_size) && length(accstem_size) == 1L &&
     accstem_size != ""){
    args[["SLEN"]] <- accstem_size
  }
  #
  args
}

#' @rdname import.tRNAdb.ce
#' @export
import.tRNAdb.ce.pattern <- function(upstream,
                                     accstem,
                                     accstem_rev,
                                     np89,
                                     dstem,
                                     dloop,
                                     dstem_rev,
                                     np26,
                                     acsteam,
                                     acloop,
                                     acstem_rev,
                                     vloop,
                                     tstem,
                                     tloop,
                                     tstem_rev,
                                     np7376,
                                     downstream,
                                     accstem_size = "both",
                                     types = c("Bacteria"),
                                     dbURL = TRNA_DB_CE_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  if(!.is_non_empty_character(accstem_size)){
    stop("'accstem_size' must be 'both', '7', or '8'")
  }
  accstem_size <- match.arg(accstem_size,
                            c("both","7","8"))
  types <- match.arg(types,
                     c(TRNA_DB_CE_TYPE,NA_character_), 
                     several.ok = TRUE)
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_ce_pattern(upstream = upstream,
                                                accstem = accstem,
                                                accstem_rev = accstem_rev,
                                                np89 = np89,
                                                dstem = dstem,
                                                dloop = dloop,
                                                dstem_rev = dstem_rev,
                                                np26 = np26,
                                                acsteam = acsteam,
                                                acloop = acloop,
                                                acstem_rev = acstem_rev,
                                                vloop = vloop,
                                                tstem = tstem,
                                                tloop = tloop,
                                                tstem_rev = tstem_rev,
                                                np7376 = np7376,
                                                downstream = downstream,
                                                accstem_size = accstem_size,
                                                types = types)
  args[["submit"]] <- "GO"
  args[["STYPE"]] <- "pattern"
  args[["VTYPE"]] <- 1L
  args <- args[c("TSEQ_1",
                 "TSEQ_2",
                 "TSEQ_3",
                 "TSEQ_4",
                 "TSEQ_5",
                 "TSEQ_6",
                 "TSEQ_7",
                 "TSEQ_8",
                 "TSEQ_9",
                 "TSEQ_10",
                 "TSEQ_11",
                 "TSEQ_12",
                 "TSEQ_13",
                 "TSEQ_14",
                 "TSEQ_15",
                 "TSEQ_16",
                 "TSEQ_17",
                 "SLEN",
                 "SEARCH_JOIN",
                 "ADV",
                 "STYPE",
                 "VTYPE",
                 "submit")]
  # get result and return as GRanges
  res <- .get_trna_db_ce_result(args = args,
                                dbURL = dbURL,
                                dbFunction = "whole_seq_search_adv.cgi",
                                verbose = verbose)
  return(.convert_tRNAdb_ce_result_to_GRanges(res))
}


# helper function --------------------------------------------------------------

# convert results from the tRNA db CE into a GRanges object
.convert_tRNAdb_ce_result_to_GRanges <- function(df){
  # convert to DataFrame is not already present
  df <- S4Vectors::DataFrame(df)
  # if empty result return empty GRanges
  if(nrow(df) == 0L){
    return(GenomicRanges::GRanges())
  }
  # 
  names(df$tRNA_seq) <- df$tRNAdb_ID
  names(df$tRNA_str) <- df$tRNAdb_ID
  # construct a valid StringSet object
  df <- .sanitize_sequences(df)
  gr <- GenomicRanges::GRanges(
    seqnames = df$tRNAdb_ID,
    ranges = IRanges::IRanges(start = rep(1,nrow(df)),
                              end = width(df$tRNA_seq)),
    strand = "*",
    df)
  gr
}

# extracting information from search list --------------------------------------

#' @importFrom httr2 request req_method req_body_form req_error req_perform 
#'    url_build
.get_trna_ce_db <- function(url, body = list(), verbose){
  req <- httr2::request(httr2::url_build(url))
  req <- httr2::req_method(req, "POST")
  req <- do.call(httr2::req_body_form, c(list(req),body))
  req <- httr2::req_error(req)
  res <- try(do.call(httr2::req_perform, 
                     list(req, verbosity = .get_verbosity(verbose))), 
             silent = TRUE)
  if(is(res,"try-error")){
    if(verbose){
      warning(res, call. = FALSE)
    } else {
      warning("tRNAdb-CE Server seems to be not available.", call. = FALSE)
    }
    return(httr2::response(status_code = 503L, 
                           url = httr2::url_build(url), 
                           method = "POST",
                           headers = c("Content-Type: text/html")))
  }
  res
}

# get main result and establish session
.get_trna_db_ce_list <- function(dbURL, dbFunction, args, verbose){
  dbURL$path <- paste0(dbURL$path, "cgi-bin/trnadb/",  dbFunction)
  .get_trna_ce_db(url = dbURL,
                  body = args,
                  verbose = verbose)
}

.get_trna_db_ce_list_page <- function(dbURL, i, verbose){
  dbURL$path <- paste0(dbURL$path, "cgi-bin/trnadb/Result")
  .get_trna_ce_db(url = dbURL,
                  body = list(position = i),
                  verbose = verbose)
}

.get_trna_db_ce_list_sequences <- function(dbURL, verbose){
  dbURL$path <- paste0(dbURL$path, "cgi-bin/trnadb/Tools")
  .get_trna_ce_db(url = dbURL,
                  body = list(),
                  verbose = verbose)
}

#' @importFrom httr2 url_parse
.get_trna_db_ce_result_detailpage <- function(dbURL, id, verbose){
  dbURL$path <- paste0(dbURL$path, "cgi-bin/trnadb/Result")
  dbURL$query <- list(ID = id)
  .get_trna_ce_db(dbURL,
                  verbose = verbose)
}

#' @importFrom httr2 resp_body_html url_parse resp_is_error
#' @importFrom IRanges CharacterList
.get_trna_db_ce_result <- function(args,
                                   dbURL,
                                   dbFunction = c("whole_seq_list_adv.cgi",
                                                  "whole_seq_search.cgi",
                                                  "whole_seq_search_adv.cgi"),
                                   verbose){
  dbURL <- httr2::url_parse(dbURL)
  dbFunction <- match.arg(dbFunction)
  # get main result and establish session
  res <- .get_trna_db_ce_list(dbURL, dbFunction, args, verbose)
  if(httr2::resp_is_error(res)){
    return(S4Vectors::DataFrame())
  }
  
  # 
  # # check result length
  # pageNumbers <- .extract_page_numbers(httr2::resp_body_html(res))
  # if(length(pageNumbers) == 0){
  #   stop("No results found.",
  #        call. = FALSE)
  # }
  # # get all pages of results
  # df <- lapply(pageNumbers,
  #              function(i){
  #                page <- .get_trna_db_list_page(dbURL, i, verbose)
  #                .extract_data_frame_from_xml_per_page(httr2::resp_body_html(page))
  #              })
  # df <- do.call(rbind,df)
  # df <- df[order(df$tRNAdb_ID),]
  # # get sequence and structure information
  # sequences <- .get_trna_db_list_sequences(dbURL, verbose)
  # sequences <- .extract_tRNAdb_sequences(httr2::resp_body_html(sequences), df)
  # sequences <- sequences[order(sequences$tRNAdb_ID),]
  # # get detail pages
  # detailpages <- lapply(df$tRNAdb_ID,
  #                       function(id){
  #                         details <- 
  #                           .get_trna_db_result_detailpage(dbURL, id, verbose)
  #                         .extract_tRNAdb_details_information(id, httr2::resp_body_html(details))
  #                       })
  # detailpages <- do.call(rbind, detailpages)
  # # make sure results match
  # if(!all(sequences$tRNAdb_ID == df$tRNAdb_ID)){
  #   stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
  #        "\nReason: unmatching tRNAdb ids.",
  #        call. = FALSE)
  # }
  # if(!all(sequences$tRNA_type == df$tRNA_type)){
  #   stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
  #        "\nReason: unmatching tRNA types.",
  #        call. = FALSE)
  # }
  # if(!all(detailpages$tRNAdb_ID == df$tRNAdb_ID)){
  #   stop("Details page returned unmatching list of tRNAdb entries.",
  #        "\nReason: unmatching tRNA ids.",
  #        call. = FALSE)
  # }
  # merge results
  df$tRNAdb_organism <- sequences$tRNAdb_organism
  df$tRNAdb_strain <- sequences$tRNAdb_strain
  df$tRNAdb_taxonomyID <- sequences$tRNAdb_taxonomyID
  df$tRNA_anticodon <- sequences$tRNA_anticodon
  df$tRNA_seq <- sequences$tRNA_seq
  df$tRNA_str <- sequences$tRNA_str
  df$tRNA_CCA.end <- sequences$tRNA_CCA.end
  df$tRNAdb_reference <- IRanges::CharacterList(as.list(detailpages$reference))
  df$tRNAdb_pmid <- IRanges::CharacterList(as.list(detailpages$pmid))
  # add additional info
  df$no <- seq_len(nrow(df))
  df$tRNA_length <- nchar(df$tRNA_seq)
  # order columns
  colOrder <- c("no",
                "tRNA_length",
                "tRNA_type",
                "tRNA_anticodon",
                "tRNA_seq",
                "tRNA_str",
                "tRNA_CCA.end",
                "tRNAdb",
                "tRNAdb_ID",
                "tRNAdb_organism",
                "tRNAdb_strain",
                "tRNAdb_taxonomyID",
                "tRNAdb_verified",
                "tRNAdb_reference",
                "tRNAdb_pmid")
  df <- df[,c(colOrder,
              colnames(df)[!(colnames(df) %in% colOrder)])]
  # save metadata
  df <- S4Vectors::DataFrame(df)
  if(dbFunction == "Blast"){
    S4Vectors::metadata(df)$BLAST <- blastRes
  }
  df
}

# html parsing -----------------------------------------------------------------
#' 
#' #' @importFrom xml2 xml_attr xml_find_all
#' .extract_ce_page_numbers <- function(xml){
#'   ans <- unique(xml2::xml_attr(xml2::xml_find_all(
#'     xml,
#'     './/td[@class="querynavtd"]//select[@name="position"]//option'),
#'     "value"))
#'   ans
#' }
#' 
#' #' @importFrom xml2 xml_attr xml_find_all
#' .extract_ce_data_frame_from_xml_per_page <- function(xml){
#'   ids <- xml2::xml_attr(xml2::xml_find_all(
#'     xml,
#'     './/tr[@class="listtabletd"]//input[@type="checkbox"][@name="selection"]'),
#'     "value")
#'   dbType <- xml2::xml_attr(xml2::xml_find_all(
#'     xml,
#'     './/tr[@class="listtabletd"]//td[2]'),
#'     "class")
#'   if(any(dbType %in% c("GRAY","gray"))){
#'     dbType[dbType %in% c("GRAY","gray")] <- "MT"
#'   }
#'   aminoacid <- as.character(xml2::xml_find_all(
#'     xml,
#'     './/tr[@class="listtabletd"]//span[@class="aminoacid"]/text()'))
#'   organism <- as.character(xml2::xml_find_all(
#'     xml,
#'     './/tr[@class="listtabletd"]//td[3]//span[@class="middle"]//a/text()'))
#'   strain <- vapply(stringr::str_split(organism," "),
#'                    function(s){
#'                      if(length(s) < 3){
#'                        return("")
#'                      }
#'                      paste(s[3:length(s)],collapse = " ")
#'                    },
#'                    character(1))
#'   strain[stringr::str_detect(strain,"\\.\\.\\.")] <- ""
#'   organism <- lapply(stringr::str_split(organism,","),"[",1)
#'   organism <- vapply(stringr::str_split(organism," "),
#'                      function(s){
#'                        paste(s[c(1,2)],collapse = " ")
#'                      },
#'                      character(1))
#'   verified <- xml2::xml_attr(xml2::xml_find_all(
#'     xml,
#'     './/tr[@class="listtabletd"]//td[2]//img'),
#'     "title")
#'   verified <- unname(TRNA_DB_VERIFIED[match(verified, names(TRNA_DB_VERIFIED))])
#'   ans <- DataFrame(tRNAdb_ID = ids,
#'                    tRNAdb = toupper(dbType),
#'                    tRNA_type = aminoacid,
#'                    tRNAdb_organism = unlist(organism),
#'                    tRNAdb_strain = unlist(strain),
#'                    tRNAdb_verified = verified)
#'   ans
#' }

#' #' @importFrom xml2 xml_text xml_find_all
#' .extract_tRNAdb_ce_details_information <- function(id, xml){
#'   reference <- ""
#'   pmid <- ""
#'   keys <- xml2::xml_text(xml2::xml_find_all(xml, '//table[@class="entrytable"][1]//div'))
#'   refkey <- which(grepl("Reference",keys))
#'   pmidkey <- which(grepl("PubMed ID",keys))
#'   if(length(refkey)){
#'     reference <- trimws(keys[refkey+1L])
#'   }
#'   if(length(pmidkey)){
#'     pmid <- trimws(keys[pmidkey+1L])
#'   }
#'   ans <- DataFrame(tRNAdb_ID = id,
#'                    reference = reference,
#'                    pmid = pmid)
#'   ans
#' }

# extract trna db information --------------------------------------------------

# 
# .extract_tRNAdb_sequences <- function(input,
#                                         df){
#   input <- stringr::str_split(input,"\n")[[1]]
#   input <- split(input[seq_len(length(input)-1)],
#                  rep(c(1,2,3),(length(input)-1)/3))
#   input[[1]] <- strsplit(input[[1]],"\\|")
#   ids <- gsub(">","",vapply(input[[1]],"[",character(1),1))
#   aminoacid <- vapply(input[[1]],"[",character(1),4)
#   anticodon <- vapply(input[[1]],"[",character(1),5)
#   taxonomyID <- vapply(input[[1]],"[",character(1),3)
#   organism <- df$tRNAdb_organism
#   strain <- df$tRNAdb_strain
#   seq <- input[[2]]
#   str <- gsub("\\(",">",gsub("\\)","<",input[[3]]))
#   if(any(stringr::str_detect(seq,"_") & df$tRNA_db != "RNA")){
#     warning("Unknown character \"_\" detected in the tRNA sequences. They ",
#             "will be removed.",
#             call. = FALSE)
#     pos <- stringr::str_locate_all(seq,"_")
#     f <- vapply(pos, function(p){nrow(p) > 0},logical(1))
#     seq[f] <- unlist(mapply(
#       function(s,p){
#         p2 <- c(0,as.numeric(as.character(p)),nchar(s)+1)
#         p <- split(p2, rep(seq_len(nrow(p)+1),nrow(p)))
#         paste(vapply(p,
#                      function(z){
#                        substr(s,z[1]+1,z[2]-1)
#                      },
#                      character(1)), collapse = "")
#       },
#       seq[f],
#       pos[f],
#       SIMPLIFY = FALSE))
#     str[f] <- unlist(mapply(
#       function(s,p){
#         p2 <- c(0,as.numeric(as.character(p)),nchar(s)+1)
#         p <- split(p2, rep(seq_len(nrow(p)+1),nrow(p)))
#         paste(vapply(p,
#                      function(z){
#                        substr(s,z[1]+1,z[2]-1)
#                      },
#                      character(1)), collapse = "")
#       },
#       str[f],
#       pos[f],
#       SIMPLIFY = FALSE))
#   }
#   # sanity check for matching sequence and structure length
#   if(!all(vapply(seq,nchar,double(1)) == vapply(str,nchar,double(1)))){
#     stop("Sequence and structure length do not match for some case:\n",
#          "\nSequences:\n",
#          paste(seq[vapply(seq,nchar,double(1)) != vapply(str,nchar,double(1))],
#                collapse = "\n"),
#          "\nStructures:\n",
#          paste(str[vapply(seq,nchar,double(1)) != vapply(str,nchar,double(1))],
#                collapse = "\n"))
#   }
#   # since apparently not all dot bracket annotations are valid, we have to 
#   # catch them. .sanitize_structures removes invalid structures. the
#   # result is now valid.
#   str <- .sanitize_structures(ids,str) # not it is a DotBracketStringSet
#   cca <- .has_CCA_end(str)
#   # create result as data.frame
#   ans <- DataFrame(tRNAdb_ID = ids,
#                    tRNA_type = aminoacid,
#                    tRNA_anticodon = anticodon,
#                    tRNAdb_organism = organism,
#                    tRNAdb_strain = strain,
#                    tRNAdb_taxonomyID = taxonomyID,
#                    tRNA_seq = seq,
#                    tRNA_str = str,
#                    tRNA_CCA.end = cca)
#   ans
# }
