### Unit tests for validation output parsing (R/validation.R) ###

context("Validation output parsing")

test_that("parse_validation_output extracts numeric statistics", {
  mock <- c(
    "Validation summary",
    "  Bias (b0):              0.0123",
    "  Dispersion (b1):        0.987",
    "  Accuracy (correlation): 0.6521",
    "  Predictive ability:     0.58",
    "done")
  s <- parse_validation_output(mock)

  # Values, not just lines — and labels like b0/b1 must not be picked up.
  expect_equal(s$bias, 0.0123)
  expect_equal(s$dispersion, 0.987)
  expect_equal(s$accuracy, 0.6521)
  expect_equal(s$predictive_ability, 0.58)

  # Matched lines are kept for inspection.
  expect_true(!is.null(s$bias_lines))
  expect_identical(s$raw, mock)
})

test_that("parse_validation_output handles missing statistics gracefully", {
  s <- parse_validation_output(c("no statistics here", "just text"))
  expect_null(s$bias)
  expect_null(s$accuracy)
  expect_identical(s$raw, c("no statistics here", "just text"))
})

test_that("parse_validation_output parses signed and exponential values", {
  s <- parse_validation_output("  bias: -1.2e-3")
  expect_equal(s$bias, -1.2e-3)
})


context("validationf90 input checks")

test_that("validationf90 rejects solutions files without a recognised header", {
  # validationf90 identifies the solutions layout from its header line; a
  # headerless file makes it abort in Fortran, so the wrapper checks first.
  d <- file.path(tempdir(), "valhdrchk")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  par_f <- file.path(d, "fake.par")
  writeLines(c("DATAFILE", "renf90.dat"), par_f)

  headerless <- file.path(d, "sol_bad")
  writeLines(c("   1   1   1   0.5", "   1   1   2   0.7"), headerless)

  expect_error(
    validationf90(par_file = par_f,
                  solutions_whole = headerless,
                  solutions_partial = headerless,
                  validation_ids = c(1L, 2L),
                  effect = 2L, dir = d),
    "solutions header recognised by validationf90")
})
