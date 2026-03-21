## Tests for heterogeneous residual variance helpers (R/hetres.R)

context("Heterogeneous residual variance helpers")

test_that("hetres_options produces class-based options", {
  result <- hetres_options(group_col = 5, n_groups = 10)
  expect_match(result[1], "^hetres_int 5 10$")
})

test_that("hetres_options produces covariate-based options", {
  result <- hetres_options(covariate_cols = c(8, 9),
                            initial = c(4.0, 0.1, 0.1))
  expect_length(result, 2)
  expect_match(result[1], "^hetres_pos 8 9$")
  expect_match(result[2], "^hetres_pol 4 0\\.1 0\\.1$")
})

test_that("hetres_options includes var_file for class-based", {
  result <- hetres_options(group_col = 5, n_groups = 3,
                            var_file = "my_hetres")
  expect_length(result, 2)
  expect_match(result[2], "^hetres_var my_hetres$")
})

test_that("hetres_options works with single covariate", {
  result <- hetres_options(covariate_cols = 10,
                            initial = c(2.5, 0.01))
  expect_match(result[1], "^hetres_pos 10$")
  expect_match(result[2], "^hetres_pol 2\\.5 0\\.01$")
})

test_that("hetres_options covariate-based without initial values", {
  result <- hetres_options(covariate_cols = 8)
  expect_length(result, 1)
  expect_match(result[1], "^hetres_pos 8$")
})

test_that("hetres_options errors with no arguments", {
  expect_error(hetres_options(), "Specify either")
})

test_that("hetres_options errors with both group and covariate", {
  expect_error(
    hetres_options(group_col = 5, n_groups = 3, covariate_cols = 8),
    "not both"
  )
})

test_that("hetres_options errors when n_groups missing for class-based", {
  expect_error(
    hetres_options(group_col = 5),
    "n_groups.*required"
  )
})
