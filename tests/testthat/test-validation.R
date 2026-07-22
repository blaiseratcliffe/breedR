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
