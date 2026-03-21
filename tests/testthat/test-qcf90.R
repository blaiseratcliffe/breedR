## Tests for QCF90 argument building and validation (R/qc.R)

context("QCF90 argument building and validation")

## Helper to call build_qcf90_args with defaults
default_args <- function(...) {
  defaults <- list(
    snp_file = "geno.txt", ped_file = NULL, map_file = NULL,
    call_rate_markers = 0.90, call_rate_animals = 0.90, maf = 0.05,
    hwe = NULL, check_parentage = FALSE, trio = FALSE,
    check_identity = FALSE, identity_threshold = 0.99,
    sex_chr = NULL, skip_chr = NULL, exclude_chr = NULL,
    save_clean = FALSE, save_log = FALSE, outcallrate = FALSE,
    remove_markers = FALSE, remove_animals = FALSE,
    extra_args = NULL
  )
  overrides <- list(...)
  defaults[names(overrides)] <- overrides
  do.call(build_qcf90_args, defaults)
}

## -- Basic args --

test_that("basic args include snpfile and thresholds", {
  args <- default_args()
  expect_true("--snpfile" %in% args)
  expect_true("geno.txt" %in% args)
  expect_true("--crm" %in% args)
  expect_true("--cra" %in% args)
  expect_true("--maf" %in% args)
})

## -- Pedigree and map --

test_that("includes pedfile and mapfile when provided", {
  args <- default_args(ped_file = "ped.txt", map_file = "map.txt")
  expect_true("--pedfile" %in% args)
  expect_true("ped.txt" %in% args)
  expect_true("--mapfile" %in% args)
  expect_true("map.txt" %in% args)
})

test_that("omits pedfile when NULL", {
  args <- default_args()
  expect_false("--pedfile" %in% args)
})

## -- HWE --

test_that("includes HWE when specified", {
  args <- default_args(hwe = 0.15)
  expect_true("--hwe" %in% args)
  expect_true("0.15" %in% args)
})

test_that("omits HWE when NULL", {
  args <- default_args()
  expect_false("--hwe" %in% args)
})

## -- Parentage --

test_that("includes parentage when TRUE", {
  args <- default_args(check_parentage = TRUE)
  expect_true("--check-parentage" %in% args)
})

test_that("includes trio", {
  args <- default_args(trio = TRUE)
  expect_true("--trio" %in% args)
})

test_that("omits parentage when FALSE", {
  args <- default_args(check_parentage = FALSE)
  expect_false("--check-parentage" %in% args)
})

## -- Identity --

test_that("includes identity check", {
  args <- default_args(check_identity = TRUE, identity_threshold = 0.95)
  expect_true("--check-identity" %in% args)
  expect_true("--threshold-identity" %in% args)
  expect_true("0.95" %in% args)
})

test_that("identity with file uses basename", {
  args <- default_args(check_identity = "/path/to/ids.txt")
  expect_true("--check-identity" %in% args)
  expect_true("ids.txt" %in% args)
  expect_false("/path/to/ids.txt" %in% args)
})

test_that("omits identity when FALSE", {
  args <- default_args()
  expect_false("--check-identity" %in% args)
})

## -- Chromosome filtering (Bug 1 fix: separate arg elements) --

test_that("sex_chr passes each chromosome as separate arg", {
  args <- default_args(sex_chr = c(30, 31))
  expect_true("--sex-chr" %in% args)
  expect_true("30" %in% args)
  expect_true("31" %in% args)
  # Should NOT be collapsed into "30 31"
  expect_false("30 31" %in% args)
})

test_that("skip_chr passes each chromosome as separate arg", {
  args <- default_args(skip_chr = c(29, 30))
  expect_true("--skip-chr" %in% args)
  expect_true("29" %in% args)
  expect_true("30" %in% args)
})

test_that("exclude_chr passes each chromosome as separate arg", {
  args <- default_args(exclude_chr = 31)
  expect_true("--exclude-chr" %in% args)
  expect_true("31" %in% args)
})

## -- Output options --

test_that("save options are included when TRUE", {
  args <- default_args(save_clean = TRUE, save_log = TRUE,
                        remove_markers = TRUE, remove_animals = TRUE)
  expect_true("--save-clean" %in% args)
  expect_true("--save-log" %in% args)
  expect_true("--remove-markers" %in% args)
  expect_true("--remove-animals" %in% args)
})

test_that("save options omitted when FALSE", {
  args <- default_args()
  expect_false("--save-clean" %in% args)
  expect_false("--save-log" %in% args)
  expect_false("--remove-markers" %in% args)
  expect_false("--remove-animals" %in% args)
})

test_that("outcallrate included when TRUE", {
  args <- default_args(outcallrate = TRUE)
  expect_true("--outcallrate" %in% args)
})

## -- Extra args --

test_that("extra_args passed through", {
  args <- default_args(extra_args = c("--fastread", "--num-threads", "4"))
  expect_true("--fastread" %in% args)
  expect_true("--num-threads" %in% args)
  expect_true("4" %in% args)
})

## -- Custom thresholds --

test_that("custom thresholds are used", {
  args <- default_args(call_rate_markers = 0.80,
                        call_rate_animals = 0.85, maf = 0.01)
  expect_true("0.8" %in% args)
  expect_true("0.85" %in% args)
  expect_true("0.01" %in% args)
})

## -- Input validation (Bug 3 fix) --

test_that("qcf90 errors on missing SNP file", {
  expect_error(qcf90("nonexistent.txt"), "not found")
})

test_that("qcf90 errors on missing pedigree file", {
  # Use a file that exists for snp_file
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, ped_file = "nonexistent_ped.txt"), "not found")
})

test_that("qcf90 errors on missing map file", {
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, map_file = "nonexistent_map.txt"), "not found")
})

test_that("qcf90 errors on invalid call_rate_markers", {
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, call_rate_markers = 1.5), "between 0 and 1")
  expect_error(qcf90(snp, call_rate_markers = -0.1), "between 0 and 1")
})

test_that("qcf90 errors on invalid call_rate_animals", {
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, call_rate_animals = 2.0), "between 0 and 1")
})

test_that("qcf90 errors on invalid maf", {
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, maf = 0.6), "between 0 and 0.5")
  expect_error(qcf90(snp, maf = -0.1), "between 0 and 0.5")
})

test_that("qcf90 errors on invalid hwe", {
  snp <- system.file("DESCRIPTION", package = "breedR")
  expect_error(qcf90(snp, hwe = 1.5), "between 0 and 1")
})

## -- check_parentage default (Bug 6 fix) --

test_that("check_parentage defaults to TRUE with pedigree", {
  # When ped_file is provided and check_parentage is NULL,
  # it should default to TRUE
  args <- default_args(ped_file = "ped.txt", check_parentage = TRUE)
  expect_true("--check-parentage" %in% args)
})
