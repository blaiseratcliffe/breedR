## Tests for expand_variance() (R/renumf90_convenience.R)

context("Variance matrix expansion")

test_that("NULL uses default diagonal", {
  v <- expand_variance(NULL, 2, c(6, 4))
  expect_equal(dim(v), c(2, 2))
  expect_equal(v[1, 1], 6)
  expect_equal(v[2, 2], 4)
  expect_equal(v[1, 2], 0)
})

test_that("scalar expands to diagonal for 1 trait", {
  v <- expand_variance(5, 1, 3)
  expect_equal(dim(v), c(1, 1))
  expect_equal(v[1, 1], 5)
})

test_that("scalar expands to diagonal for 2 traits", {
  v <- expand_variance(5, 2, c(3, 4))
  expect_equal(dim(v), c(2, 2))
  expect_equal(v[1, 1], 5)
  expect_equal(v[2, 2], 5)
  expect_equal(v[1, 2], 0)
})

test_that("scalar expands to diagonal for 3 traits", {
  v <- expand_variance(10, 3, c(1, 2, 3))
  expect_equal(dim(v), c(3, 3))
  expect_equal(diag(v), c(10, 10, 10))
  expect_equal(v[1, 2], 0)
})

test_that("vector expands to diagonal", {
  v <- expand_variance(c(5, 3), 2, c(1, 1))
  expect_equal(dim(v), c(2, 2))
  expect_equal(v[1, 1], 5)
  expect_equal(v[2, 2], 3)
  expect_equal(v[1, 2], 0)
})

test_that("vector of length 3 for 3 traits", {
  v <- expand_variance(c(10, 8, 6), 3, c(1, 1, 1))
  expect_equal(diag(v), c(10, 8, 6))
})

test_that("matrix passes through unchanged", {
  m <- matrix(c(5, 2, 2, 3), 2, 2)
  v <- expand_variance(m, 2, c(1, 1))
  expect_identical(v, m)
})

test_that("3x3 matrix passes through", {
  m <- matrix(c(10, 3, 2, 3, 8, 1, 2, 1, 6), 3, 3)
  v <- expand_variance(m, 3, c(1, 1, 1))
  expect_identical(v, m)
})

test_that("wrong vector length errors", {
  expect_error(expand_variance(c(1, 2, 3), 2, c(1, 1)),
               "scalar.*vector.*matrix")
})

test_that("wrong matrix dimensions error", {
  m <- matrix(1:9, 3, 3)
  expect_error(expand_variance(m, 2, c(1, 1)),
               "scalar.*vector.*matrix")
})

test_that("1x1 matrix for 1 trait passes through", {
  m <- matrix(7)
  v <- expand_variance(m, 1, 3)
  expect_equal(v[1, 1], 7)
})
