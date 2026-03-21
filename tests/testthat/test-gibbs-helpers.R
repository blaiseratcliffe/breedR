## Tests for Gibbs sampling option helpers (R/gibbs.R)

context("Gibbs sampling option helpers")

## -- cat --

test_that("build_gibbs_options generates cat option", {
  result <- build_gibbs_options(cat = c(0, 2))
  expect_equal(result, "cat 0 2")
})

test_that("build_gibbs_options handles multi-category", {
  result <- build_gibbs_options(cat = c(0, 0, 2, 5))
  expect_equal(result, "cat 0 0 2 5")
})

test_that("build_gibbs_options rejects negative cat", {
  expect_error(build_gibbs_options(cat = c(0, -1)), "non-negative")
})

## -- thresholds --

test_that("build_gibbs_options generates thresholds", {
  result <- build_gibbs_options(thresholds = c(0.0, 1.0, 2.0))
  expect_match(result, "^thresholds 0 1 2$")
})

test_that("build_gibbs_options rejects non-numeric thresholds", {
  expect_error(build_gibbs_options(thresholds = "abc"), "numeric")
})

## -- prior --

test_that("build_gibbs_options generates prior", {
  result <- build_gibbs_options(prior = c(5, 2, -1, 5))
  expect_equal(result, "prior 5 2 -1 5")
})

## -- seed --

test_that("build_gibbs_options generates seed", {
  result <- build_gibbs_options(seed = c(123, -432))
  expect_equal(result, "seed 123 -432")
})

test_that("build_gibbs_options rejects wrong length seed", {
  expect_error(build_gibbs_options(seed = 123), "length 2")
})

## -- fixed_var --

test_that("build_gibbs_options generates fixed_var mean", {
  result <- build_gibbs_options(fixed_var = "mean")
  expect_equal(result, "fixed_var mean")
})

test_that("build_gibbs_options generates fixed_var all", {
  result <- build_gibbs_options(fixed_var = "all")
  expect_equal(result, "fixed_var all")
})

test_that("build_gibbs_options fixed_var with specific effects", {
  result <- build_gibbs_options(fixed_var = "mean", save_effects = c(1, 2, 3))
  expect_equal(result, "fixed_var mean 1 2 3")
})

test_that("build_gibbs_options rejects invalid fixed_var", {
  expect_error(build_gibbs_options(fixed_var = "bad"), "mean.*or.*all")
})

## -- save_samples (solution) --

test_that("build_gibbs_options generates solution mean", {
  result <- build_gibbs_options(save_samples = "mean")
  expect_equal(result, "solution mean")
})

test_that("build_gibbs_options generates solution all with effects", {
  result <- build_gibbs_options(save_samples = "all", save_effects = c(2, 3))
  expect_equal(result, "solution all 2 3")
})

test_that("build_gibbs_options rejects invalid save_samples", {
  expect_error(build_gibbs_options(save_samples = "bad"), "mean.*or.*all")
})

## -- cont --

test_that("build_gibbs_options generates cont", {
  result <- build_gibbs_options(cont = TRUE)
  expect_equal(result, "cont 1")
})

test_that("build_gibbs_options skips cont when FALSE", {
  result <- build_gibbs_options(cont = FALSE)
  expect_length(result, 0)
})

## -- save_halfway --

test_that("build_gibbs_options generates save_halfway_samples", {
  result <- build_gibbs_options(save_halfway = 5000)
  expect_equal(result, "save_halfway_samples 5000")
})

test_that("build_gibbs_options rejects non-positive save_halfway", {
  expect_error(build_gibbs_options(save_halfway = -1), "positive")
})

## -- hetres_int --

test_that("build_gibbs_options generates hetres_int", {
  result <- build_gibbs_options(hetres_int = list(col = 5, n = 10))
  expect_equal(result, "hetres_int 5 10")
})

test_that("build_gibbs_options rejects invalid hetres_int", {
  expect_error(build_gibbs_options(hetres_int = list(col = 5)),
               "col.*and.*n")
})

## -- residual --

test_that("build_gibbs_options generates residual", {
  result <- build_gibbs_options(residual_var = 1)
  expect_equal(result, "residual 1")
})

## -- combined --

test_that("build_gibbs_options generates multiple options", {
  result <- build_gibbs_options(
    cat = c(0, 2),
    seed = c(42, -99),
    save_samples = "mean",
    save_halfway = 1000
  )
  expect_length(result, 4)
  expect_match(result[1], "^cat 0 2$")
  expect_match(result[2], "^seed 42 -99$")
  expect_match(result[3], "^solution mean$")
  expect_match(result[4], "^save_halfway_samples 1000$")
})

test_that("build_gibbs_options returns empty for no args", {
  result <- build_gibbs_options()
  expect_length(result, 0)
})
