## Tests for SNP file I/O (R/genomic.R: write_snp_file, read_snp_file)

context("SNP file I/O")

## -- write_snp_file() integer format --

test_that("write_snp_file writes correct integer format", {
  geno <- matrix(c(0, 1, 2, 1, 0, 2), nrow = 2, byrow = TRUE)
  ids <- c("A1", "A2")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f)
  lines <- readLines(f)

  expect_length(lines, 2)
  expect_match(lines[1], "A1.*012")
  expect_match(lines[2], "A2.*102")
})

test_that("write_snp_file produces fixed-width columns", {
  geno <- matrix(c(0, 1, 2, 0, 1, 2, 0, 1, 2), nrow = 3)
  ids <- c("1", "999999", "42")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f)
  lines <- readLines(f)

  # All lines should have the same length
  expect_equal(length(unique(nchar(lines))), 1)
})

test_that("write_snp_file replaces NA with missing code", {
  geno <- matrix(c(0, NA, 2, 1), nrow = 2)
  ids <- c("A1", "A2")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f)
  lines <- readLines(f)

  # NA should become 5 (default missing)
  expect_match(lines[1], "02")
  expect_match(lines[2], "51")
})

test_that("write_snp_file custom missing value", {
  geno <- matrix(c(0, NA, 2, 1), nrow = 2)
  ids <- c("A1", "A2")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f, missing_value = 9L)
  lines <- readLines(f)
  expect_match(lines[2], "91")
})

test_that("write_snp_file errors on dimension mismatch", {
  geno <- matrix(1:6, nrow = 2)
  expect_error(write_snp_file(geno, c("A1"), tempfile()),
               "must match")
})

## -- write_snp_file() fractional format --

test_that("write_snp_file writes correct fractional format", {
  geno <- matrix(c(0.5, 1.12, 0.25, 1.50), nrow = 2)
  ids <- c("A1", "A2")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f, fractional = TRUE)
  lines <- readLines(f)

  expect_match(lines[1], "0\\.500\\.25")
  expect_match(lines[2], "1\\.121\\.50")
})

test_that("fractional format has fixed-width columns", {
  geno <- matrix(runif(12), nrow = 3)
  ids <- c("1", "99999", "42")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f, fractional = TRUE)
  lines <- readLines(f)
  expect_equal(length(unique(nchar(lines))), 1)
})

## -- read_snp_file() integer format --

test_that("read_snp_file roundtrips integer genotypes", {
  geno_orig <- matrix(c(0, 1, 2, 5, 1, 0, 2, 1, 0), nrow = 3)
  ids_orig <- c("animal1", "animal2", "animal3")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno_orig, ids_orig, f)
  result <- read_snp_file(f)

  expect_equal(result$ids, ids_orig)
  expect_equal(result$geno, geno_orig)
})

test_that("read_snp_file roundtrips with numeric IDs", {
  geno_orig <- matrix(c(0, 2, 1, 1), nrow = 2)
  ids_orig <- c("80", "8014")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno_orig, ids_orig, f)
  result <- read_snp_file(f)

  expect_equal(result$ids, ids_orig)
  expect_equal(result$geno, geno_orig)
})

## -- read_snp_file() fractional format --

test_that("read_snp_file roundtrips fractional genotypes", {
  geno_orig <- matrix(c(2.00, 0.50, 1.00, 1.12, 0.25, 1.50), nrow = 2)
  ids_orig <- c("A1", "A2")
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno_orig, ids_orig, f, fractional = TRUE)
  result <- read_snp_file(f, fractional = TRUE)

  expect_equal(result$ids, ids_orig)
  expect_equal(result$geno, geno_orig, tolerance = 0.005)
})

## -- read_snp_file() dimensions --

test_that("read_snp_file returns correct dimensions", {
  geno <- matrix(sample(0:2, 50, replace = TRUE), nrow = 5)
  ids <- paste0("ID", 1:5)
  f <- tempfile()
  on.exit(unlink(f))

  write_snp_file(geno, ids, f)
  result <- read_snp_file(f)

  expect_equal(nrow(result$geno), 5)
  expect_equal(ncol(result$geno), 10)
  expect_length(result$ids, 5)
})
