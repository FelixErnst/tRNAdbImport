#' @include tRNAdbImport.R
NULL

#' @name import.tRNAdb
#' @aliases import.tRNAdb import.mttRNAdb import.tRNAdb.id import.mttRNAdb.id
#' import.tRNAdb.blast import.mttRNAdb.blast tRNAdb2GFF
#' 
#' @title Importing information from the tRNA db as GRanges object
#' 
#' @description
#' title
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
#' @param input a GRanges object which passes the \code{istRNAdbGRanges} check
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
#' import.tRNAdb(organism = "Saccharomyces cerevisiae",
#'               aminoacids = c("Phe","Ala"))
#' import.tRNAdb.id(tdbID = "tdbD00000785")
#' import.tRNAdb.blast(blastSeq =
#' "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCA")
#' import.mttRNAdb(organism = "Bos taurus",
#'                 aminoacids = c("Phe","Ala"))
#' import.mttRNAdb.id(mtdbID = "mtdbD00000900")
NULL

#' @rdname import.tRNAdb
#' @export
TRNA_DB_URL <- "http://trna.bioinf.uni-leipzig.de/"
#' @rdname import.tRNAdb
#' @export
TRNA_DB_URL_MT <- "http://mttrna.bioinf.uni-leipzig.de/mt"


.assemble_args_for_tRNA_db_search <- 
function(organism, strain, taxonomyID, aminoacids, anticodons, sequences,
       structures, reference, comment, pubmed, genes, tdbID, blastSequence,
       database, origin){
  # base values
  args <- list(search = 0)
  #
  if(!missing(origin) && length(origin) == 1L && !is.na(origin) && origin != ""){
    args[TRNA_DB_ORIGIN[names(TRNA_DB_ORIGIN) %in% origin]] <- "on"
  }
  #
  if(!missing(database) && length(database) == 1L && !is.na(database) && 
     database != ""){
    args[["database"]] <- database
  }
  #
  if(!missing(organism) && length(organism) == 1L && !is.na(organism) && 
     organism != ""){
    args[["org"]] <- organism
  }
  #
  if(!missing(strain) && length(strain) == 1L && !is.na(strain) && 
     strain != ""){
    args[["strain"]] <- strain
  }
  #
  if(!missing(taxonomyID) && length(taxonomyID) == 1L && !is.na(taxonomyID) && 
     taxonomyID != ""){
    args[["TAX_ID"]] <- taxonomyID
  }
  #
  if(!missing(aminoacids) && all(!is.na(aminoacids)) && all(aminoacids != "")){
    aminoacids <- as.list(aminoacids)
    names(aminoacids) <- rep("aminoacid",length(aminoacids))
    args <- append(args,
                   aminoacids)
  }
  #
  if(!missing(anticodons) && all(!is.na(anticodons)) && all(anticodons != "")){
    args[["antis"]] <- paste(anticodons,collapse = ",")
  }
  #
  if(!missing(reference) && length(reference) == 1L && !is.na(reference) && 
     reference != ""){
    args[["ref"]] <- reference
  }
  #
  if(!missing(comment) && length(comment) == 1L && !is.na(comment) && 
     comment != ""){
    args[["comment"]] <- comment
  }
  #
  if(!missing(pubmed) && all(!is.na(pubmed)) && all(pubmed != "")){
    args[["pubmed"]] <- pubmed
  }
  #
  if(!missing(genes) && all(!is.na(genes)) && all(genes != "")){
    args[["gen"]] <- genes
  }
  #
  if(!missing(tdbID) && all(!is.na(tdbID)) && all(tdbID != "")){
    args[["searchID"]] <- tdbID
  }
  #
  if(!missing(blastSequence) && all(!is.na(blastSequence)) && 
     all(blastSequence != "")){
    args[["sequence"]] <- blastSequence
  }
  #
  args
}

#' @rdname import.tRNAdb
#' @export
import.tRNAdb.id <- function(tdbID, database = c("DNA", "RNA"),
                             origin = c("allothers", "plastid", "mitochondrial"),
                             dbURL = TRNA_DB_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  if(!.is_non_empty_character(tdbID)){
    stop("'tdbID' must contain only non empty character values.")
  }
  database <- match.arg(database[1], c(TRNA_DB_TYPE,NA_character_))
  origin <- match.arg(origin[1], c(names(TRNA_DB_ORIGIN),NA_character_))
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_search(tdbID = tdbID,
                                            database = database,
                                            origin = origin)
  args[["submit"]] <- "search ID"
  # get result and return as GRanges
  res <- .get_trna_db_result(args = args, dbURL = dbURL, dbFunction = "Search",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}

#' @rdname import.tRNAdb
#' @export
import.mttRNAdb.id <- function(mtdbID, dbURL = TRNA_DB_URL_MT,
                               verbose = FALSE){
  import.tRNAdb.id(tdbID = mtdbID,
                   database = NA_character_,
                   origin = NA_character_,
                   dbURL = dbURL,
                   verbose = verbose)
}

#' @rdname import.tRNAdb
#' @export
import.tRNAdb.blast <- function(blastSeq, database = c("DNA",  "RNA"),
                                origin = c("allothers", "plastid", "mitochondrial"),
                                dbURL = TRNA_DB_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  blastSeq <- as.character(blastSeq) # in case it is a single DNAString*
  if(!.is_non_empty_string(blastSeq)){
    stop("'blastSeq' must be a single non empty character value.")
  }
  database <- match.arg(database[1], c(TRNA_DB_TYPE,NA_character_))
  origin <- match.arg(origin[1], c(names(TRNA_DB_ORIGIN),NA_character_))
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_search(blastSequence = blastSeq,
                                            database = database,
                                            origin = origin)
  args[["submit"]] <- "BLAST"
  args <- args[c("database",
                 "search",
                 "sequence",
                 "submit")]
  # get result and return as GRanges
  res <- .get_trna_db_result(args = args, dbURL = dbURL, dbFunction = "Blast",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}

#' @rdname import.tRNAdb
#' @export
import.tRNAdb <- function(organism = "", strain = "", taxonomyID = "",
                          aminoacids = "", anticodons = "",  sequences = list(),
                          structures = list(), reference = "",  comment  = "",
                          pubmed = "", genes = "",  database = c("DNA", "RNA"),
                          origin = c("allothers", "plastid", "mitochondrial"),
                          dbURL = TRNA_DB_URL, verbose = FALSE){
  # input check
  if(!.is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE")
  }
  .checkValueValidity(
    aminoacids, 
    c(Biostrings::AMINO_ACID_CODE[names(Biostrings::AMINO_ACID_CODE) %in% 
                                    unique(Biostrings::GENETIC_CODE)],
      "Pyr",
      "Sec",
      "Ini",
      ""))
  .checkValueValidity(anticodons, c(names(Biostrings::GENETIC_CODE),
                                    ""))
  database <- match.arg(database[1], c(TRNA_DB_TYPE,NA_character_))
  origin <- match.arg(origin[1], c(names(TRNA_DB_ORIGIN),NA_character_))
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_search(organism = organism,
                                            strain = strain,
                                            taxonomyID = taxonomyID,
                                            aminoacids = aminoacids,
                                            anticodons = anticodons,
                                            sequences = sequences,
                                            structures = structures,
                                            reference = reference,
                                            comment  = comment,
                                            pubmed = pubmed,
                                            genes = genes,
                                            database = database,
                                            origin = origin)
  args[["submit"]] <- "search database"
  # get result and return as GRanges
  res <- .get_trna_db_result(args = args, dbURL = dbURL, dbFunction = "Search",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}

#' @rdname import.tRNAdb
#' @export
import.mttRNAdb <- function(organism = "",  strain = "",  taxonomyID = "",
                            aminoacids = "", anticodons = "", sequences = list(),
                            structures = list(), reference = "", comment = "",
                            pubmed = "", genes = "", dbURL = TRNA_DB_URL_MT,
                            verbose = FALSE){
  import.tRNAdb(organism = organism,
                strain = strain,
                taxonomyID = taxonomyID,
                aminoacids = aminoacids,
                anticodons = anticodons,
                sequences = sequences,
                structures = structures,
                reference = reference,
                comment = comment,
                pubmed = pubmed,
                genes = genes,
                database = NA_character_,
                origin = NA_character_,
                dbURL = dbURL,
                verbose = verbose)
}

# helper function --------------------------------------------------------------

# convert results from the tRNA db into a GRanges object
.convert_tRNAdb_result_to_GRanges <- function(df){
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


.handle_trna_db_connection_error <- function(res){
  
  
}

.get_verbose <- function(verbose){
  if(verbose){
    2L
  } else {
    0L
  }
}

#' @importFrom httr2 request req_method req_body_form req_error req_perform 
#'    url_build
.get_trna_db <- function(url, body = list(), verbose){
  req <- httr2::request(httr2::url_build(url))
  req <- httr2::req_method(req, "POST")
  req <- do.call(httr2::req_body_form, c(list(req),body))
  req <- httr2::req_error(req)
  res <- try(do.call(httr2::req_perform, 
                     list(req, verbosity = .get_verbose(verbose))), 
             silent = TRUE)
  if(is(res,"try-error")){
    if(verbose){
      warning(res, call. = FALSE)
    } else {
      warning("tRNAdb Server seems to be not available.", call. = FALSE)
    }
    return(httr2::response(status_code = 503L, 
                           url = httr2::url_build(url), 
                           method = "POST",
                           headers = c("Content-Type: text/html")))
  }
  res
}

# get main result and establish session
.get_trna_db_list <- function(dbURL, dbFunction, args, verbose){
  dbURL$path <- paste0(dbURL$path, "DataOutput/",  dbFunction)
  .get_trna_db(url = dbURL,
               body = args,
               verbose = verbose)
}

.get_trna_db_list_page <- function(dbURL, i, verbose){
  dbURL$path <- paste0(dbURL$path, "DataOutput/Result")
  .get_trna_db(url = dbURL,
               body = list(position = i),
               verbose = verbose)
}

.get_trna_db_list_sequences <- function(dbURL, verbose){
  dbURL$path <- paste0(dbURL$path, "DataOutput/Tools")
  .get_trna_db(url = dbURL,
               body = list("allnormalchecked" = "on",
                           "function" = "fastastruct",
                           "sel" = ""),
               verbose = verbose)
}

#' @importFrom httr2 resp_body_html
#' @importFrom xml2 xml_text xml_find_all
.get_trna_db_list_blast <- function(res, dbFunction, verbose){
  blastRes <- ""
  if(dbFunction == "Blast"){
    blastRes <- xml2::xml_text(xml2::xml_find_all(
      httr2::resp_body_html(res), './/tt'))
    if(verbose){
      message(blastRes)
    }
  }
  blastRes
}

#' @importFrom httr2 url_parse
.get_trna_db_result_detailpage <- function(dbURL, id, verbose){
  dbURL$path <- paste0(dbURL$path, "DataOutput/Result")
  dbURL$query <- list(ID = id)
  .get_trna_db(dbURL,
               verbose = verbose)
}

#' @importFrom httr2 resp_body_html url_parse resp_is_error
#' @importFrom IRanges CharacterList
.get_trna_db_result <- function(args, dbURL, dbFunction = c("Search","Blast"),
                                verbose){
  dbURL <- httr2::url_parse(dbURL)
  dbFunction <- match.arg(dbFunction)
  # get main result and establish session
  res <- .get_trna_db_list(dbURL, dbFunction, args, verbose)
  if(httr2::resp_is_error(res)){
    return(S4Vectors::DataFrame())
  }
  # output blast results if BLAST search
  blastRes <- .get_trna_db_list_blast(res, dbFunction, verbose)
  # check result length
  pageNumbers <- .extract_page_numbers(httr2::resp_body_html(res))
  if(length(pageNumbers) == 0){
    stop("No results found.",
         call. = FALSE)
  }
  # get all pages of results
  df <- lapply(pageNumbers,
               function(i){
                 page <- .get_trna_db_list_page(dbURL, i, verbose)
                 .extract_data_frame_from_xml_per_page(httr2::resp_body_html(page))
               })
  df <- do.call(rbind,df)
  df <- df[order(df$tRNAdb_ID),]
  # get sequence and structure information
  sequences <- .get_trna_db_list_sequences(dbURL, verbose)
  sequences <- .extract_tRNAdb_sequences(httr2::resp_body_html(sequences), df)
  sequences <- sequences[order(sequences$tRNAdb_ID),]
  # get detail pages
  detailpages <- lapply(df$tRNAdb_ID,
                        function(id){
                          details <- 
                            .get_trna_db_result_detailpage(dbURL, id, verbose)
                          .extract_tRNAdb_details_information(id, httr2::resp_body_html(details))
                        })
  detailpages <- do.call(rbind, detailpages)
  # make sure results match
  if(!all(sequences$tRNAdb_ID == df$tRNAdb_ID)){
    stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
         "\nReason: unmatching tRNAdb ids.",
         call. = FALSE)
  }
  if(!all(sequences$tRNA_type == df$tRNA_type)){
    stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
         "\nReason: unmatching tRNA types.",
         call. = FALSE)
  }
  if(!all(detailpages$tRNAdb_ID == df$tRNAdb_ID)){
    stop("Details page returned unmatching list of tRNAdb entries.",
         "\nReason: unmatching tRNA ids.",
         call. = FALSE)
  }
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

#' @importFrom xml2 xml_attr xml_find_all
.extract_page_numbers <- function(xml){
  ans <- unique(xml2::xml_attr(xml2::xml_find_all(
    xml,
    './/td[@class="querynavtd"]//select[@name="position"]//option'),
    "value"))
  ans
}

#' @importFrom xml2 xml_attr xml_find_all
.extract_data_frame_from_xml_per_page <- function(xml){
  ids <- xml2::xml_attr(xml2::xml_find_all(
    xml,
    './/tr[@class="listtabletd"]//input[@type="checkbox"][@name="selection"]'),
    "value")
  dbType <- xml2::xml_attr(xml2::xml_find_all(
    xml,
    './/tr[@class="listtabletd"]//td[2]'),
    "class")
  if(any(dbType %in% c("GRAY","gray"))){
    dbType[dbType %in% c("GRAY","gray")] <- "MT"
  }
  aminoacid <- as.character(xml2::xml_find_all(
    xml,
    './/tr[@class="listtabletd"]//span[@class="aminoacid"]/text()'))
  organism <- as.character(xml2::xml_find_all(
    xml,
    './/tr[@class="listtabletd"]//td[3]//span[@class="middle"]//a/text()'))
  strain <- vapply(stringr::str_split(organism," "),
                   function(s){
                     if(length(s) < 3){
                       return("")
                     }
                     paste(s[3:length(s)],collapse = " ")
                   },
                   character(1))
  strain[stringr::str_detect(strain,"\\.\\.\\.")] <- ""
  organism <- lapply(stringr::str_split(organism,","),"[",1)
  organism <- vapply(stringr::str_split(organism," "),
                     function(s){
                       paste(s[c(1,2)],collapse = " ")
                     },
                     character(1))
  verified <- xml2::xml_attr(xml2::xml_find_all(
    xml,
    './/tr[@class="listtabletd"]//td[2]//img'),
    "title")
  verified <- unname(TRNA_DB_VERIFIED[match(verified, names(TRNA_DB_VERIFIED))])
  ans <- DataFrame(tRNAdb_ID = ids,
                   tRNAdb = toupper(dbType),
                   tRNA_type = aminoacid,
                   tRNAdb_organism = unlist(organism),
                   tRNAdb_strain = unlist(strain),
                   tRNAdb_verified = verified)
  ans
}

#' @importFrom xml2 xml_text xml_find_all
.extract_tRNAdb_details_information <- function(id, xml){
  reference <- ""
  pmid <- ""
  keys <- xml2::xml_text(xml2::xml_find_all(xml, '//table[@class="entrytable"][1]//div'))
  refkey <- which(grepl("Reference",keys))
  pmidkey <- which(grepl("PubMed ID",keys))
  if(length(refkey)){
    reference <- trimws(keys[refkey+1L])
  }
  if(length(pmidkey)){
    pmid <- trimws(keys[pmidkey+1L])
  }
  ans <- DataFrame(tRNAdb_ID = id,
                   reference = reference,
                   pmid = pmid)
  ans
}

# formating sequences ----------------------------------------------------------

.pos_letters <- function(x,chrs){
  lapply(chrs,
         function(chr){
           stringr::str_locate_all(x, chr)
         })
}
.sanitize_structures <- function(ids,str){
  open <- .pos_letters(str,Structstrings::STRUCTURE_OPEN_CHR)
  close <- .pos_letters(str,Structstrings::STRUCTURE_CLOSE_CHR)
  lengthOpen <- lapply(open,function(z){lengths(z)})
  lengthClose <- lapply(close,function(z){lengths(z)})
  lengthMatch <- lapply(
    seq_along(lengthOpen),
    function(i){
      which(unlist(lengthOpen[[i]]) != unlist(lengthClose[[i]]))
    })
  f <- unique(unlist(lengthMatch))
  if(length(f) > 0L){
    warning("Result from the tRNAdb contain invalid dot bracket annotation.\n",
            "The following tRNAdb ids contain the invalid structure ",
            "information: '",
            paste(ids[f],
                  collapse = "', '"),
            "'")
    strLengths <- unlist(lapply(str[f], width))
    newStr <- unlist(lapply(strLengths,
                            function(len){
                              paste(rep(".",len),collapse = "")
                            }))
    str[f] <- newStr
  }
  # this checks for validity
  str <- Structstrings::DotBracketStringSet(str)
}
.sanitize_sequences <- function(df){
  seqs <- df$tRNA_seq
  f_dna <- which(df$tRNAdb == "DNA")
  f_rna <- which(df$tRNAdb == "RNA" | df$tRNAdb == "MT")
  seq_dna <- seqs[f_dna]
  seq_rna <- seqs[f_rna]
  if(length(seq_dna) > 0){
    seq_dna <- Biostrings::DNAStringSet(seq_dna)
  } else {
    seq_dna <- NULL
  }
  if(length(seq_rna) > 0){
    seq_rna <- Modstrings::sanitizeFromtRNAdb(seq_rna)
    seq_rna <- gsub("_","",seq_rna) # removes the insertion character
    seq_modrna <- Modstrings::ModRNAStringSet(seq_rna)
    seq_rna <- as(seq_modrna,"RNAStringSet")
    seq_dna_test <- as(seq_rna,"DNAStringSet")
    # if the ModRNAstringSet is actually a DNAStringset
    if(all(as.character(seq_dna_test) == as.character(seq_modrna))){
      seq_rna <- seq_dna_test
    } else {
      # if the ModRNAstringSet does contain modifications keep the ModRNAStringSet
      if(!all(as.character(seq_rna) == as.character(seq_modrna))){
        seq_rna <- seq_modrna
      }
    }
    rm(seq_dna_test)
    rm(seq_modrna)
  } else {
    seq_rna <- NULL
  }
  if(is(seq_rna,"ModRNAStringSet") || is(seq_rna,"RNAStringSet")){
    seq_dna <- as(seq_dna,"RNAStringSet")
  }
  if(is(seq_rna,"ModRNAStringSet")){
    seq_dna <- as(seq_dna,"ModRNAStringSet")
  }
  if(!is.null(seq_dna) &
     !is.null(seq_rna) &
     class(seq_rna) != class(seq_dna)){
    stop("Something went wrong.")
  }
  seqs <- list(seq_dna,seq_rna) 
  seqs <- do.call(c,
                  seqs[!vapply(seqs,is.null,logical(1))])
  seqs <- seqs[c(f_dna,f_rna)]
  df$tRNA_seq <- seqs
  df
}


# extract trna db information --------------------------------------------------

.has_CCA_end <- function(structures){
  strList <- getBasePairing(structures)
  vapply(strList,
         function(str){
           end <- max(str$pos)
           # the last three nucleotides must be unpaired
           all(str[str$pos %in% (end-2):end,]$forward == 0)
         },
         logical(1))
}

.extract_tRNAdb_sequences <- function(input,
                                        df){
  input <- stringr::str_split(input,"\n")[[1]]
  input <- split(input[seq_len(length(input)-1)],
                 rep(c(1,2,3),(length(input)-1)/3))
  input[[1]] <- strsplit(input[[1]],"\\|")
  ids <- gsub(">","",vapply(input[[1]],"[",character(1),1))
  aminoacid <- vapply(input[[1]],"[",character(1),4)
  anticodon <- vapply(input[[1]],"[",character(1),5)
  taxonomyID <- vapply(input[[1]],"[",character(1),3)
  organism <- df$tRNAdb_organism
  strain <- df$tRNAdb_strain
  seq <- input[[2]]
  str <- gsub("\\(",">",gsub("\\)","<",input[[3]]))
  if(any(stringr::str_detect(seq,"_") & df$tRNA_db != "RNA")){
    warning("Unknown character \"_\" detected in the tRNA sequences. They ",
            "will be removed.",
            call. = FALSE)
    pos <- stringr::str_locate_all(seq,"_")
    f <- vapply(pos, function(p){nrow(p) > 0},logical(1))
    seq[f] <- unlist(mapply(
      function(s,p){
        p2 <- c(0,as.numeric(as.character(p)),nchar(s)+1)
        p <- split(p2, rep(seq_len(nrow(p)+1),nrow(p)))
        paste(vapply(p,
                     function(z){
                       substr(s,z[1]+1,z[2]-1)
                     },
                     character(1)), collapse = "")
      },
      seq[f],
      pos[f],
      SIMPLIFY = FALSE))
    str[f] <- unlist(mapply(
      function(s,p){
        p2 <- c(0,as.numeric(as.character(p)),nchar(s)+1)
        p <- split(p2, rep(seq_len(nrow(p)+1),nrow(p)))
        paste(vapply(p,
                     function(z){
                       substr(s,z[1]+1,z[2]-1)
                     },
                     character(1)), collapse = "")
      },
      str[f],
      pos[f],
      SIMPLIFY = FALSE))
  }
  # sanity check for matching sequence and structure length
  if(!all(vapply(seq,nchar,double(1)) == vapply(str,nchar,double(1)))){
    stop("Sequence and structure length do not match for some case:\n",
         "\nSequences:\n",
         paste(seq[vapply(seq,nchar,double(1)) != vapply(str,nchar,double(1))],
               collapse = "\n"),
         "\nStructures:\n",
         paste(str[vapply(seq,nchar,double(1)) != vapply(str,nchar,double(1))],
               collapse = "\n"))
  }
  # since apparently not all dot bracket annotations are valid, we have to 
  # catch them. .sanitize_structures removes invalid structures. the
  # result is now valid.
  str <- .sanitize_structures(ids,str) # not it is a DotBracketStringSet
  cca <- .has_CCA_end(str)
  # create result as data.frame
  ans <- DataFrame(tRNAdb_ID = ids,
                   tRNA_type = aminoacid,
                   tRNA_anticodon = anticodon,
                   tRNAdb_organism = organism,
                   tRNAdb_strain = strain,
                   tRNAdb_taxonomyID = taxonomyID,
                   tRNA_seq = seq,
                   tRNA_str = str,
                   tRNA_CCA.end = cca)
  ans
}

# convert result to GFF format compatible result

#' @rdname import.tRNAdb
#' @export
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
