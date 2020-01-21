library(tRNAdbImport)

context("tRNAdbImport")
test_that("tRNAdbImport:",{
  httptest::skip_if_disconnected(url = "http://trna.bioinf.uni-leipzig.de/")
  gr <- import.tRNAdb(organism = "Saccharomyces cerevisiae",
                      aminoacids = c("Phe","Ala"))
  expect_equal(c("no","tRNA_length","tRNA_type","tRNA_anticodon","tRNA_seq",
                 "tRNA_str","tRNA_CCA.end","tRNAdb","tRNAdb_ID",
                 "tRNAdb_organism","tRNAdb_strain","tRNAdb_taxonomyID",
                 "tRNAdb_verified"),colnames(mcols(gr)))
  #
  expect_true(istRNAdbGRanges(gr))
  #
  length <- as.numeric(S4Vectors::mcols(gr)$tRNA_length)
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_seq))
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_str))
  #
  expect_type(mcols(gr)$no, "integer")
  expect_type(mcols(gr)$tRNA_length, "integer")
  expect_type(mcols(gr)$tRNA_type, "character")
  expect_type(mcols(gr)$tRNA_anticodon, "character")
  expect_type(mcols(gr)$tRNA_seq, "S4")
  expect_type(mcols(gr)$tRNA_str, "S4")
  expect_type(mcols(gr)$tRNA_CCA.end, "logical")
  expect_type(mcols(gr)$tRNAdb_ID, "character")
  expect_type(mcols(gr)$tRNAdb_organism, "character")
  expect_type(mcols(gr)$tRNAdb_strain, "character")
  expect_type(mcols(gr)$tRNAdb_taxonomyID, "character")
  expect_type(mcols(gr)$tRNAdb_verified, "logical")
  #
  gff <- tRNAdb2GFF(gr)
  expect_type(mcols(gff)$source, "integer")
  expect_true(is.factor(mcols(gff)$source))
  expect_type(mcols(gff)$type, "integer")
  expect_true(is.factor(mcols(gff)$type))
  expect_type(mcols(gff)$score, "integer")
  expect_type(mcols(gff)$phase, "integer")
  expect_type(mcols(gff)$ID, "character")
  # 
  expect_type(mcols(gff)$no, "integer")
  expect_type(mcols(gff)$tRNA_length, "integer")
  expect_type(mcols(gff)$tRNA_type, "character")
  expect_type(mcols(gff)$tRNA_anticodon, "character")
  expect_type(mcols(gff)$tRNA_seq, "character")
  expect_type(mcols(gff)$tRNA_str, "character")
  expect_type(mcols(gff)$tRNA_CCA.end, "logical")
  expect_type(mcols(gff)$tRNAdb_ID, "character")
  expect_type(mcols(gff)$tRNAdb_organism, "character")
  expect_type(mcols(gff)$tRNAdb_strain, "character")
  expect_type(mcols(gff)$tRNAdb_taxonomyID, "character")
  expect_type(mcols(gff)$tRNAdb_verified, "logical")
  #
  mcols(gr)$tRNA_seq <- NULL
  expect_warning(expect_false(istRNAdbGRanges(gr)),
                 "Input GRanges object does not meet the requirements of the")
  expect_warning(expect_false(.check_trnadb_granges(data.frame())),
                 "Input is not a GRanges object")
  # other general input tests
  gr <- import.tRNAdb.id(tdbID = "tdbD00000785")
  expect_true(istRNAdbGRanges(gr))
  gr <- import.tRNAdb.blast(blastSeq =
    "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCA")
  expect_true(istRNAdbGRanges(gr))
  expect_equal(names(S4Vectors::metadata(mcols(gr))),"BLAST")
  gr <- import.mttRNAdb(organism = "Bos taurus", aminoacids = c("Phe","Ala"))
  expect_true(istRNAdbGRanges(gr))
  gr <- import.mttRNAdb.id(mtdbID = "mtdbD00000900")
  expect_true(istRNAdbGRanges(gr))
  gr <- import.tRNAdb(organism = "Saccharomyces cerevisiae", anticodons = c("GAA","CAT"))
  expect_true(istRNAdbGRanges(gr))
  # verbose output
  expect_output(import.tRNAdb.id(tdbID = "tdbD00000785", verbose = TRUE))
  expect_output(import.tRNAdb.blast(
    blastSeq = "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCA",
    verbose = TRUE))
})

context("input failure and warning tests")
test_that("input failure and warning test:",{
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search()
  expect_type(actual,"list")
  expect_named(actual,"search")
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(list(test = "12"))
  expect_type(actual,"list")
  expect_named(actual,c("search","org"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(strain = "test")
  expect_named(actual,c("search","strain"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(taxonomyID = "test")
  expect_named(actual,c("search","TAX_ID"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(anticodons = "test")
  expect_named(actual,c("search","antis"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(reference = "test")
  expect_named(actual,c("search","ref"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(comment = "test")
  expect_named(actual,c("search","comment"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(pubmed = "test")
  expect_named(actual,c("search","pubmed"))
  actual <- tRNAdbImport:::.assemble_args_for_tRNA_db_search(genes = "comment")
  expect_named(actual,c("search","gen"))
})

context("input failure and warning tests")
test_that("input failure and warning test:",{
  expect_error(
    tRNAdbImport:::.has_CCA_end(),
    'argument "structures" is missing'
  )
  expect_error(
    tRNAdbImport:::.has_CCA_end("")
  )
  expect_error(import.tRNAdb.id(tdbID = "tdbD0000078500"))
  expect_error(tRNAdbImport:::.checkValueValidity("a", c("b","c")),
               "'a' must be one of the following values: 'b', 'c'")
  
})
