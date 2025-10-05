
# httr2 helper functions

.get_verbosity <- function(verbose){
  if(verbose){
    2L
  } else {
    0L
  }
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
