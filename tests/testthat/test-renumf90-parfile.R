## Tests for RENUMF90 parameter file builder (R/renumf90.R)

context("RENUMF90 parameter file builder")

test_that("build_renumf90_parfile generates basic single-trait file", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = 3,
    residual_variance = 10,
    effects = list(
      list(pos = 1, type = "cross", form = "alpha")
    ),
    fields_passed = NULL, weights = NULL, snp_file = NULL,
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  expect_true("DATAFILE" %in% result)
  expect_true("data.txt" %in% result)
  expect_true("TRAITS" %in% result)
  expect_true("3" %in% result)
  expect_true(any(grepl("1 cross alpha", result)))
  expect_true("RESIDUAL_VARIANCE" %in% result)
  expect_true("10" %in% result)
})

test_that("build_renumf90_parfile handles multi-trait", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = c(5, 6),
    residual_variance = matrix(c(10, 3, 3, 8), 2, 2),
    effects = list(
      list(pos = c(2, 2), type = "cross", form = "alpha")
    ),
    fields_passed = NULL, weights = NULL, snp_file = NULL,
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  expect_true(any(grepl("^5 6$", result)))
  expect_true(any(grepl("2 2 cross alpha", result)))
  expect_true(any(grepl("^10 3$", result)))
  expect_true(any(grepl("^3 8$", result)))
})

test_that("build_renumf90_parfile handles different effects per trait", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = c(5, 6),
    residual_variance = matrix(c(10, 3, 3, 8), 2, 2),
    effects = list(
      list(pos = c(2, 2), type = "cross", form = "alpha"),
      list(pos = c(3, 0), type = "cov"),
      list(pos = c(0, 4), type = "cross", form = "alpha")
    ),
    fields_passed = NULL, weights = NULL, snp_file = NULL,
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  expect_true(any(grepl("3 0 cov", result)))
  expect_true(any(grepl("0 4 cross alpha", result)))
})

test_that("build_renumf90_parfile includes animal random effect", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = 3,
    residual_variance = 10,
    effects = list(
      list(pos = 1, type = "cross", form = "alpha",
           random = "animal",
           file = "pedigree.txt",
           covariances = 5)
    ),
    fields_passed = NULL, weights = NULL, snp_file = NULL,
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  expect_true("RANDOM" %in% result)
  expect_true("animal" %in% result)
  expect_true("FILE" %in% result)
  expect_true("pedigree.txt" %in% result)
  expect_true("(CO)VARIANCES" %in% result)
  expect_true("5" %in% result)
})

test_that("build_renumf90_parfile includes optional keywords", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = 3,
    residual_variance = 10,
    effects = list(
      list(pos = 1, type = "cross", form = "alpha",
           random = "animal",
           file = "ped.txt",
           file_pos = c(1, 2, 3, 0, 4),
           covariances = 5,
           optional = c("pe", "mat"),
           covariances_pe = 2)
    ),
    fields_passed = c(2, 4), weights = 7,
    snp_file = NULL, ped_depth = 0, inbreeding = "pedigree",
    upg_type = NULL, gen_int = NULL, skip_header = 2,
    options = c("alpha_size 30")
  )

  expect_true(any(grepl("^2 4$", result)))
  expect_true(any(grepl("^7$", result)))
  expect_true("SKIP_HEADER" %in% result)
  expect_true("2" %in% result)
  expect_true("FILE_POS" %in% result)
  expect_true(any(grepl("1 2 3 0 4", result)))
  expect_true("PED_DEPTH" %in% result)
  expect_true("INBREEDING" %in% result)
  expect_true("pedigree" %in% result)
  expect_true("OPTIONAL" %in% result)
  expect_true(any(grepl("pe mat", result)))
  expect_true(any(grepl("OPTION alpha_size 30", result)))
})

test_that("build_renumf90_parfile includes SNP_FILE for animal effects", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = 3,
    residual_variance = 10,
    effects = list(
      list(pos = 1, type = "cross", form = "alpha",
           random = "animal",
           file = "ped.txt",
           covariances = 5)
    ),
    fields_passed = NULL, weights = NULL,
    snp_file = "genotypes.txt",
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  expect_true("SNP_FILE" %in% result)
  expect_true("genotypes.txt" %in% result)
})

test_that("build_renumf90_parfile default residual when NULL", {
  result <- build_renumf90_parfile(
    datafile = "data.txt",
    traits = 3,
    residual_variance = NULL,
    effects = list(
      list(pos = 1, type = "cross", form = "alpha")
    ),
    fields_passed = NULL, weights = NULL, snp_file = NULL,
    ped_depth = NULL, inbreeding = NULL, upg_type = NULL,
    gen_int = NULL, skip_header = NULL, options = NULL
  )

  # Should have a default residual variance
  res_idx <- which(result == "RESIDUAL_VARIANCE")
  expect_true(length(res_idx) == 1)
  expect_equal(result[res_idx + 1], "1")
})
