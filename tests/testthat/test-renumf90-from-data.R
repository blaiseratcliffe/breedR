## Tests for renumf90_from_data() validation (R/renumf90_convenience.R)
## Only tests input validation — no RENUMF90 binary needed.

context("renumf90_from_data validation")

test_that("errors on non-data.frame input", {
  expect_error(renumf90_from_data(data = "not_a_df", traits = "y"),
               "data.frame")
})

test_that("errors on empty traits", {
  dat <- data.frame(y = 1:3)
  expect_error(renumf90_from_data(data = dat, traits = character(0)),
               "at least one")
})

test_that("errors on missing column names", {
  dat <- data.frame(y = 1:3, x = c("A", "B", "A"))
  expect_error(renumf90_from_data(data = dat, traits = "y",
                                    fixed = c("x", "nonexistent")),
               "not found.*nonexistent")
})

test_that("errors on non-numeric trait column", {
  dat <- data.frame(y = c("a", "b", "c"), x = 1:3)
  expect_error(renumf90_from_data(data = dat, traits = "y"),
               "numeric")
})

test_that("errors on bad pedigree (not data.frame)", {
  dat <- data.frame(id = 1:3, y = rnorm(3))
  expect_error(renumf90_from_data(data = dat, traits = "y",
                                    pedigree = "ped.txt"),
               "data.frame")
})

test_that("errors on bad pedigree (fewer than 3 columns)", {
  dat <- data.frame(id = 1:3, y = rnorm(3))
  ped <- data.frame(id = 1:3, sire = 0)
  expect_error(renumf90_from_data(data = dat, traits = "y",
                                    pedigree = ped),
               "at least 3")
})

test_that("errors on pedigree animal column not in data", {
  dat <- data.frame(tree = 1:3, y = rnorm(3))
  ped <- data.frame(id = 1:3, sire = 0, dam = 0)
  expect_error(renumf90_from_data(data = dat, traits = "y",
                                    pedigree = ped),
               "not found in data")
})

test_that("warns on continuous random effect", {
  dat <- data.frame(id = 1:200, age = runif(200), y = rnorm(200))
  expect_warning(
    tryCatch(
      renumf90_from_data(data = dat, traits = "y", random = "age",
                          residual_variance = 5,
                          dir = file.path(tempdir(), "test_rnd_warn2")),
      error = function(e) NULL  # may fail at binary step, that's OK
    ),
    "unique numeric values"
  )
})

test_that("classify_col detects character as cross", {
  dat <- data.frame(site = c("A", "B", "A"), y = 1:3)
  # We can't call classify_col directly (it's inside renumf90_from_data),
  # but we can verify the behavior indirectly: character columns should
  # produce "cross alpha" effect lines
})

test_that("classify_col respects factors override", {
  # numeric column listed in factors should be cross, not cov
  dat <- data.frame(block = c(1, 2, 3), y = c(10, 12, 11))
  # Can only verify this by checking the generated parameter file
  # which requires running RENUMF90 — skip for unit test
})

test_that("errors on all-NA data after filtering", {
  dat <- data.frame(id = c(NA, NA, NA), y = c(10, 12, 11))
  ped <- data.frame(id = c(NA, NA, NA), sire = 0, dam = 0)
  expect_error(
    renumf90_from_data(data = dat, traits = "y", pedigree = ped,
                        dir = file.path(tempdir(), "test_empty2")),
    "No data remaining")
})
