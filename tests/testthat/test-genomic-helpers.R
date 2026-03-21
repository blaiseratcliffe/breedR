## Tests for genomic helper functions (R/genomic.R)

context("Genomic helper functions")

## -- check_genomic() --

test_that("check_genomic requires snp_file", {
  expect_error(check_genomic(list()), "snp_file.*required")
})

test_that("check_genomic rejects non-existent snp_file", {
  expect_error(check_genomic(list(snp_file = "nonexistent.txt")),
               "not found")
})

test_that("check_genomic rejects non-existent map_file", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), map_file = "no.txt")),
    "not found")
})

test_that("check_genomic sets correct defaults", {
  result <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR")))
  expect_equal(result$whichG, 1L)
  expect_equal(result$tunedG, 2L)
  expect_equal(result$AlphaBeta, c(0.95, 0.05))
  expect_equal(result$minfreq, 0.05)
  expect_equal(result$callrate, 0.90)
  expect_equal(result$callrateAnim, 0.90)
  expect_equal(result$verify_parentage, 3L)
  expect_false(result$saveG)
  expect_false(result$saveA22)
})

test_that("check_genomic rejects invalid whichG", {
  expect_error(check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), whichG = 5)),
               "whichG.*must be 1, 2, or 3")
})

test_that("check_genomic rejects invalid tunedG", {
  expect_error(check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), tunedG = 7)),
               "tunedG")
})

test_that("check_genomic rejects invalid AlphaBeta", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), AlphaBeta = 0.95)),
    "length 2")
})

test_that("check_genomic rejects invalid callrate", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), callrate = 1.5)),
    "between 0 and 1")
})

test_that("check_genomic rejects invalid verify_parentage", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), verify_parentage = 5)),
    "verify_parentage")
})

test_that("check_genomic accepts valid extra_options", {
  result <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                extra_options = c("no_quality_control")))
  expect_equal(result$extra_options, "no_quality_control")
})

test_that("check_genomic rejects non-character extra_options", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), extra_options = 42)),
    "character vector")
})

## -- build_genomic_options() --

test_that("build_genomic_options generates correct options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR")))
  opts <- build_genomic_options(genomic)

  expect_true(any(grepl("^SNP_file DESCRIPTION$", opts)))
  expect_true(any(grepl("^whichG 1$", opts)))
  expect_true(any(grepl("^tunedG 2$", opts)))
  expect_true(any(grepl("^AlphaBeta 0\\.95 0\\.05$", opts)))
  expect_true(any(grepl("^minfreq 0\\.05$", opts)))
  expect_true(any(grepl("^callrate 0\\.9$", opts)))
  expect_true(any(grepl("^verify_parentage 3$", opts)))
})

test_that("build_genomic_options includes map_file when provided", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 map_file = system.file("NAMESPACE", package = "breedR")))
  opts <- build_genomic_options(genomic)
  expect_true(any(grepl("^map_file", opts)))
})

test_that("build_genomic_options includes save options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 saveG = TRUE, saveA22 = TRUE))
  opts <- build_genomic_options(genomic)
  expect_true("saveG" %in% opts)
  expect_true("saveA22" %in% opts)
})

test_that("build_genomic_options includes extra_options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 extra_options = c("no_quality_control",
                                                   "thrStopCorAG 0.0")))
  opts <- build_genomic_options(genomic)
  expect_true("no_quality_control" %in% opts)
  expect_true("thrStopCorAG 0.0" %in% opts)
})
