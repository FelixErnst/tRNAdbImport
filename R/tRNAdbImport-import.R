#' @include tRNAdbImport.R
NULL

#' @name import.tRNAdb
#' @aliases import.tRNAdb import.mttRNAdb import.tRNAdb.id import.mttRNAdb.id
#' import.tRNAdb.blast import.mttRNAdb.blast
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
#' @param verbose whether to report verbose information from the httr calls
#'
#' @return a GRanges object containing the information from the tRNA db
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors DataFrame
#' @importFrom stringr str_detect str_split str_locate_all
#' 
#' @export
#'
#' @examples
#' \donttest{
#' import.tRNAdb(organism = "Saccharomyces cerevisiae",
#'               aminoacids = c("Phe","Ala"),
#'               anticodons = c("GAA"))
#' import.tRNAdb.id(tdbID = "tdbD00000785")
#' import.tRNAdb.blast(blastSeq = 
#' "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCA")
#' import.mttRNAdb(organism = "Bos taurus")
#' }
NULL

#' @rdname import.tRNAdb
#' @export
import.tRNAdb.id <- function(tdbID,
                             database = "DNA",
                             origin = "allothers",
                             dbURL = TRNA_DB_URL,
                             verbose = FALSE){
  # input check
  assertive::assert_is_a_bool(verbose)
  assertive::assert_all_are_non_empty_character(tdbID)
  .checkValueValidity(database, TRNA_DB_TYPE)
  .checkValueValidity(origin, c(names(TRNA_DB_ORIGIN),NA))
  # assemble arguments
  args <- .assemble_args_for_tRNA_db_search(tdbID = tdbID,
                                            database = database,
                                            origin = origin)
  args[["submit"]] <- "search ID"
  # get result and return as GRanges
  res <- .get_trna_db_result(args = args,
                             dbURL = dbURL,
                             dbFunction = "Search",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}
#' @rdname import.tRNAdb
#' @export
import.mttRNAdb.id <- function(mtdbID,
                               dbURL = TRNA_DB_URL_MT,
                               verbose = FALSE){
  import.tRNAdb.id(tdbID = mtdbID,
                   database = NA,
                   origin = NA,
                   dbURL = dbURL,
                   verbose = verbose)
}

#' @rdname import.tRNAdb
#' @export
import.tRNAdb.blast <- function(blastSeq,
                                database = "DNA",
                                origin = "allothers",
                                dbURL = TRNA_DB_URL,
                                verbose = FALSE){
  # input check
  assertive::assert_is_a_bool(verbose)
  blastSeq <- as.character(blastSeq) # in case it is a single DNAString*
  assertive::assert_is_a_non_empty_string(blastSeq)
  .checkValueValidity(database, TRNA_DB_TYPE)
  .checkValueValidity(origin, c(names(TRNA_DB_ORIGIN),NA))
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
  res <- .get_trna_db_result(args = args,
                             dbURL = dbURL,
                             dbFunction = "Blast",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}
# #' @rdname import.tRNAdb
# #' @export
# import.mttRNAdb.blast <- function(blastSeq,
#                                   dbURL = TRNA_DB_URL_MT,
#                                   verbose = FALSE){
#   import.tRNAdb.blast(blastSeq = blastSeq,
#                       database = NA,
#                       origin = NA,
#                       dbURL = dbURL,
#                       verbose = verbose)
# }

#' @rdname import.tRNAdb
#' @export
import.tRNAdb <- function(organism = "",
                          strain = "",
                          taxonomyID = "",
                          aminoacids = "",
                          anticodons = "",
                          sequences = list(),
                          structures = list(),
                          reference = "",
                          comment  = "",
                          pubmed = "",
                          genes = "",
                          database = "DNA",
                          origin = "allothers",
                          dbURL = TRNA_DB_URL,
                          verbose = FALSE){
  # input check
  assertive::assert_is_a_bool(verbose)
  .checkValueValidity(aminoacids, c(Biostrings::AMINO_ACID_CODE[names(Biostrings::AMINO_ACID_CODE) %in% unique(Biostrings::GENETIC_CODE)],
                                    "Pyr",
                                    "Sec",
                                    "Ini",
                                    ""))
  .checkValueValidity(anticodons, c(names(Biostrings::GENETIC_CODE),
                                    ""))
  .checkValueValidity(database, TRNA_DB_TYPE)
  .checkValueValidity(origin, c(names(TRNA_DB_ORIGIN),NA))
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
  res <- .get_trna_db_result(args = args,
                             dbURL = dbURL,
                             dbFunction = "Search",
                             verbose = verbose)
  return(.convert_tRNAdb_result_to_GRanges(res))
}
#' @rdname import.tRNAdb
#' @export
import.mttRNAdb <- function(organism = "",
                            strain = "",
                            taxonomyID = "",
                            aminoacids = "",
                            anticodons = "",
                            sequences = list(),
                            structures = list(),
                            reference = "",
                            comment = "",
                            pubmed = "",
                            genes = "",
                            dbURL = TRNA_DB_URL_MT,
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
                database = NA,
                origin = NA,
                dbURL = dbURL,
                verbose = verbose)
}

# convert results from the tRNA db into a GRanges object
.convert_tRNAdb_result_to_GRanges <- function(df){
  df <- S4Vectors::DataFrame(df)
  gr <- GenomicRanges::GRanges(
    seqnames = df$tRNA_dbID,
    ranges = IRanges::IRanges(start = rep(1,nrow(df)),
                              end = nchar(df$tRNA_seq)),
    strand = "*",
    df)
  if(all(df$tRNA_db != "RNA")){
    gr$tRNA_seq <- Biostrings::DNAStringSet(gr$tRNA_seq)
  } else {
    warning("The result potentially contains modified DNA and RNA nucleotides.",
            "\nDownstream analysis might be limited.",
            call. = FALSE)
  }
  gr
}

.get_trna_db_result <- function(args,
                                dbURL,
                                dbFunction = "Search",
                                verbose){
  # get main result and establish session
  if(verbose){
    res <- httr::POST(paste0(httr::modify_url(dbURL),
                             "DataOutput/",
                             dbFunction),
                      body = args,
                      encode = "form",
                      httr::verbose())
  } else {
    res <- httr::POST(paste0(httr::modify_url(dbURL),
                             "DataOutput/",
                             dbFunction),
                      body = args,
                      encode = "form")
  }
  # output blast results if BLAST search
  if(dbFunction == "Blast" && verbose){
    message(xml_text(xml2::xml_find_all(content(res), './/tt')))
  }
  # check result length
  pageNumbers <- .extract_page_numbers(httr::content(res))
  if(length(pageNumbers) == 0){
    stop("No results found.",
         call. = FALSE)
  }
  # get all pages of results
  df <- lapply(pageNumbers,
               function(i){
                 if(verbose){
                   page <- httr::POST(paste0(httr::modify_url(dbURL),
                                             "DataOutput/Result"),
                                      body = list(position = i),
                                      encode = "form",
                                      httr::verbose())
                 } else {
                   page <- httr::POST(paste0(httr::modify_url(dbURL),
                                             "DataOutput/Result"),
                                      body = list(position = i),
                                      encode = "form")
                 }
                 ans <- 
                   .extract_data_frame_from_xml_per_page(httr::content(page))
                 ans
               })
  df <- do.call(rbind,df)
  df <- df[order(df$tRNA_dbID),]
  # get sequence and structure information
  if(verbose){
    seqs <- httr::POST(paste0(httr::modify_url(dbURL),
                              "DataOutput/Tools"),
                       body = list("allnormalchecked" = "on",
                                   "function" = "fastastruct",
                                   "sel" = ""),
                       encode = "form",
                       httr::verbose())
  } else {
    seqs <- httr::POST(paste0(httr::modify_url(dbURL),
                              "DataOutput/Tools"),
                       body = list("allnormalchecked" = "on",
                                   "function" = "fastastruct",
                                   "sel" = ""),
                       encode = "form")
  }
  seqs <- .extract_sequences_and_structures(httr::content(seqs),
                                            df)
  seqs <- seqs[order(seqs$tRNA_dbID),]
  # make sure results match
  if(!all(seqs$tRNA_dbID == df$tRNA_dbID)){
    stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
         "\nReason: unmatching tRNAdb ids.",
         call. = FALSE)
  }
  if(!all(seqs$tRNA_type == df$tRNA_type)){
    stop("Function 'fastastruct' returned unmatching list of tRNAdb entries.",
         "\nReason: unmatching tRNA types.",
         call. = FALSE)
  }
  # merge results
  df$tRNA_organism <- seqs$tRNA_organism
  df$tRNA_strain <- seqs$tRNA_strain
  df$tRNA_taxonomyID <- seqs$tRNA_taxonomyID
  df$tRNA_anticodon <- seqs$tRNA_anticodon
  df$tRNA_seq <- seqs$tRNA_seq
  df$tRNA_str <- seqs$tRNA_str
  df$tRNA_CCA.end <- seqs$tRNA_CCA.end
  # add additional info
  df$no <- 1:nrow(df)
  df$tRNA_length <- nchar(df$tRNA_seq)
  # order columns
  colOrder <- c("no",
                "tRNA_length",
                "tRNA_type",
                "tRNA_anticodon",
                "tRNA_seq",
                "tRNA_str",
                "tRNA_CCA.end",
                "tRNA_db",
                "tRNA_dbID",
                "tRNA_organism",
                "tRNA_strain",
                "tRNA_taxonomyID",
                "tRNA_verified")
  df <- df[,c(colOrder,
              colnames(df)[!(colnames(df) %in% colOrder)])]
  df
}

.extract_page_numbers <- function(xml){
  ans <- unique(xml2::xml_attr(xml2::xml_find_all(xml, './/td[@class="querynavtd"]//select[@name="position"]//option'),"value"))
  ans
}

.extract_data_frame_from_xml_per_page <- function(xml){
  ids <- xml2::xml_attr(xml2::xml_find_all(xml, './/tr[@class="listtabletd"]//input[@type="checkbox"][@name="selection"]'),"value")
  dbType <- xml2::xml_attr(xml2::xml_find_all(xml, './/tr[@class="listtabletd"]//td[2]'),"class")
  if(any(dbType == "GRAY")){
    dbType[dbType == "GRAY"] <- "mt"
  }
  aminoacid <- as.character(xml2::xml_find_all(xml, './/tr[@class="listtabletd"]//span[@class="aminoacid"]/text()'))
  organism <- as.character(xml2::xml_find_all(xml, './/tr[@class="listtabletd"]//td[3]//span[@class="middle"]//a/text()'))
  organism <- lapply(stringr::str_split(organism,","),"[",1)
  strain <- vapply(stringr::str_split(organism," "),
                   function(s){
                     if(length(s) < 3){
                       return("")
                     }
                     paste(s[3:length(s)],collapse = " ")
                   },
                   character(1))
  strain[stringr::str_detect(strain,"\\.\\.\\.")] <- ""
  organism <- vapply(stringr::str_split(organism," "),
                     function(s){
                       paste(s[1:2],collapse = " ")
                     },
                     character(1))
  verified <- xml2::xml_attr(xml2::xml_find_all(xml, './/tr[@class="listtabletd"]//td[2]//img'),"title")
  verified <- unname(TRNA_DB_VERIFIED[match(verified, names(TRNA_DB_VERIFIED))])
  ans <- data.frame(tRNA_dbID = ids,
                    tRNA_db = toupper(dbType),
                    tRNA_type = aminoacid,
                    tRNA_organism = unlist(organism),
                    tRNA_strain = unlist(strain),
                    tRNA_verified = verified,
                    stringsAsFactors = FALSE)
  ans
}

.extract_sequences_and_structures <- function(input,
                                              df){
  input <- stringr::str_split(input,"\n")[[1]]
  input <- split(input[1:(length(input)-1)],
                rep(c(1,2,3),(length(input)-1)/3))
  input[[1]] <- strsplit(input[[1]],"\\|")
  ids <- gsub(">","",vapply(input[[1]],"[",character(1),1))
  aminoacid <- vapply(input[[1]],"[",character(1),4)
  anticodon <- vapply(input[[1]],"[",character(1),5)
  taxonomyID <- vapply(input[[1]],"[",character(1),3)
  organism <- gsub("_"," ",vapply(input[[1]],"[",character(1),2))
  strain <- vapply(stringr::str_split(organism," "),
                   function(s){
                     if(length(s) < 3){
                       return("")
                     }
                     paste(s[3:length(s)],collapse = " ")
                   },
                   character(1))
  organism <- vapply(stringr::str_split(organism," "),
                     function(s){
                       paste(s[1:2],collapse = " ")
                     },
                     character(1))
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
        p <- split(p2, rep(1:(nrow(p)+1),nrow(p)))
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
        p <- split(p2, rep(1:(nrow(p)+1),nrow(p)))
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
  # catch them
  tryCatch({
    cca <- .has_CCA_end(str)
  },
  error = function(e){
    stop("Result from the tRNAdb contain invalid dot bracket annotation.\n",
         e,
         "\ntRNAdb ids: '",
         paste(ids[as.numeric(gsub("'","",stringr::str_extract(e,"'[0-9]++'")))],
               collapse = "', '"),
         "'",
         call. = FALSE)
  })
  # create result as data.frame
  ans <- data.frame(tRNA_dbID = ids,
                    tRNA_type = aminoacid,
                    tRNA_anticodon = anticodon,
                    tRNA_organism = organism,
                    tRNA_strain = strain,
                    tRNA_taxonomyID = taxonomyID,
                    tRNA_seq = seq,
                    tRNA_str = str,
                    tRNA_CCA.end = cca,
                    stringsAsFactors = FALSE)
  ans
}

.has_CCA_end <- function(structures){
  strList <- tRNA::getBasePairing(structures)
  vapply(strList,
         function(str){
           end <- max(str$pos)
           all(str[str$pos %in% (end-2):end,]$forward == 0)
         },
         logical(1))
}

.assemble_args_for_tRNA_db_search <- function(organism,
                                              strain,
                                              taxonomyID,
                                              aminoacids,
                                              anticodons,
                                              sequences,
                                              structures,
                                              reference,
                                              comment,
                                              pubmed,
                                              genes,
                                              tdbID,
                                              blastSequence,
                                              database,
                                              origin){
  # base values
  args <- list(search = 0)
  #
  if(!missing(origin) && !is.na(origin) && origin != ""){
    args[TRNA_DB_ORIGIN[names(TRNA_DB_ORIGIN) %in% origin]] <- "on"
  }
  #
  if(!missing(database) && !is.na(database) && database != ""){
    args[["database"]] <- database
  }
  #
  if(!missing(organism) && !is.na(organism) && organism != ""){
    args[["org"]] <- organism
  }
  #
  if(!missing(strain) && !is.na(strain) && strain != ""){
    args[["strain"]] <- strain
  }
  #
  if(!missing(taxonomyID) && !is.na(taxonomyID) && taxonomyID != ""){
    args[["TAX_ID"]] <- taxonomyID
  }
  #
  if(!missing(aminoacids) && !is.na(aminoacids) && aminoacids != ""){
    args[["aminoacid"]] <- list(aminoacids)
  }
  #
  if(!missing(anticodons) && !is.na(anticodons) && anticodons != ""){
    args[["antis"]] <- list(anticodons)
  }
  #
  if(!missing(reference) && !is.na(reference) && reference != ""){
    args[["ref"]] <- reference
  }
  #
  if(!missing(comment) && !is.na(comment) && comment != ""){
    args[["comment"]] <- comment
  }
  #
  if(!missing(pubmed) && !is.na(pubmed) && pubmed != ""){
    args[["pubmed"]] <- pubmed
  }
  #
  if(!missing(genes) && !is.na(genes) && genes != ""){
    args[["gen"]] <- genes
  }
  #
  if(!missing(tdbID) && !is.na(tdbID) && tdbID != ""){
    args[["searchID"]] <- tdbID
  }
  #
  if(!missing(blastSequence) && !is.na(blastSequence) && blastSequence != ""){
    args[["sequence"]] <- blastSequence
  }
  #
  args
}