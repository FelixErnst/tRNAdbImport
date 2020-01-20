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
  expect_warning(expect_false(istRNAdbGRanges(gr)))
  # other general input tests
  gr <- import.tRNAdb.id(tdbID = "tdbD00000785")
  expect_true(istRNAdbGRanges(gr))
  gr <- import.tRNAdb.blast(blastSeq =
    "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCA")
  expect_true(istRNAdbGRanges(gr))
  gr <- import.mttRNAdb(organism = "Bos taurus", aminoacids = c("Phe","Ala"))
  expect_true(istRNAdbGRanges(gr))
  gr <- import.mttRNAdb.id(mtdbID = "mtdbD00000900")
  expect_true(istRNAdbGRanges(gr))
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
})
