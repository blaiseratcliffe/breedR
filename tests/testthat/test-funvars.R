context("Functions of variance components")

## AIREMLF90's OPTION se_covar_function reports three numbers per function:
## the plug-in value at the REML solution ('mean'), and the mean and standard
## deviation of its Monte Carlo sampling draws. summary() must report the
## plug-in value as the Estimate, so that it agrees with the variance
## components printed directly above it (issue #13).

test_that("funvars is parsed as a 3-row matrix", {

  res <- load_res("ped_ar")

  expect_true(is.matrix(res$funvars))
  expect_identical(dim(res$funvars), c(3L, 1L))
  expect_identical(rownames(res$funvars),
                   c('mean', 'sample mean', 'sample sd'))
  expect_identical(colnames(res$funvars), "Heritability")

  ## multi-trait model: one column per trait
  res_mt <- load_res("mt")
  expect_identical(dim(res_mt$funvars), c(3L, 2L))
  expect_identical(rownames(res_mt$funvars),
                   c('mean', 'sample mean', 'sample sd'))
})


test_that("the Estimate is the plug-in value, not the Monte Carlo mean", {

  res <- load_res("ped_ar")

  ## the plug-in heritability implied by the reported variance components
  v <- res$var[, 1]
  h2 <- unname(v['genetic'] / sum(v))

  ## the 'mean' row is that value; the 'sample mean' row is not
  expect_equal(unname(res$funvars['mean', 1]), h2, tolerance = 1e-4)
  expect_false(isTRUE(all.equal(unname(res$funvars['sample mean', 1]), h2,
                                tolerance = 1e-5)))

  tbl <- funvars_table(res$funvars)
  expect_identical(dim(tbl), c(1L, 3L))
  expect_identical(colnames(tbl), c("Estimate", "Sample Mean", "S.E."))
  expect_identical(rownames(tbl), "Heritability")
  expect_equal(unname(tbl[1, "Estimate"]), h2, tolerance = 1e-4)
  expect_equal(unname(tbl[1, "Sample Mean"]),
               unname(res$funvars['sample mean', 1]))
  expect_equal(unname(tbl[1, "S.E."]), unname(res$funvars['sample sd', 1]))
})


test_that("summary() prints the plug-in estimate alongside the sampling draws", {

  res <- load_res("ped_ar")
  v <- res$var[, 1]
  h2 <- unname(v['genetic'] / sum(v))

  expect_output(print(summary(res)), "Heritability")
  expect_output(print(summary(res)), "Estimate")
  expect_output(print(summary(res)), "Sample Mean")

  out <- capture.output(print(summary(res)))
  hline <- grep("^Heritability", out, value = TRUE)
  expect_length(hline, 1L)

  nums <- as.numeric(strsplit(trimws(sub("^Heritability", "", hline)),
                              " +")[[1]])
  expect_length(nums, 3L)

  ## the first printed column is the plug-in estimate (0.28458), and
  ## specifically not the Monte Carlo sample mean (0.28415)
  expect_equal(nums[1], h2, tolerance = 1e-4)
  expect_equal(nums[2], unname(res$funvars['sample mean', 1]), tolerance = 1e-4)
  expect_equal(nums[3], unname(res$funvars['sample sd', 1]), tolerance = 1e-4)
})


test_that("divergence between the estimate and the sample mean is flagged", {

  fv <- function(...) {
    ans <- cbind(c(...))
    rownames(ans) <- c('mean', 'sample mean', 'sample sd')
    colnames(ans) <- 'Heritability'
    ans
  }

  ## well-determined variance components: the two means agree
  ## res_ped_ar, gap of 0.005 sampling SD
  expect_false(unname(funvars_divergent(fv(0.284580, 0.284150, 0.093236))))

  ## globulus (issue #13), gap of 0.020 sampling SD
  expect_false(unname(funvars_divergent(fv(0.190280, 0.188500, 0.088017))))

  ## harvey (issue #13), gap of 0.13 sampling SD on a quantity bounded in [0, 1]
  expect_true(unname(funvars_divergent(fv(0.71741, 0.62372, 0.70984))))

  ## a degenerate sampling SD must not produce NA/NaN
  expect_false(unname(funvars_divergent(fv(0.10950, 0.10950, 0))))
  expect_true(unname(funvars_divergent(fv(0.10950, 0.20950, 0))))

  ## missing values are not divergences
  expect_false(unname(funvars_divergent(fv(0.5, NA, 0.1))))

  ## no variance functions requested (EM, or no genetic effect):
  ## parse_functions() returns list()
  expect_identical(funvars_divergent(list()), logical(0))

  ## names are carried through, one per function, so the note can name them.
  ## a single column is the case that needs watching: indexing a one-column
  ## matrix by row drops the column name
  expect_identical(funvars_divergent(fv(0.71741, 0.62372, 0.70984)),
                   c(Heritability = TRUE))

  fv2 <- cbind(a = c(0.5, 0.5, 0.1), b = c(0.71741, 0.62372, 0.70984))
  rownames(fv2) <- c('mean', 'sample mean', 'sample sd')
  expect_identical(funvars_divergent(fv2), c(a = FALSE, b = TRUE))
})


test_that("summary() explains why the estimate and the sample mean diverge", {

  res <- load_res("ped_ar")

  ## no note when the variance components are well determined
  expect_false(any(grepl("Note", capture.output(print(summary(res))))))

  ## substitute the divergent harvey values from issue #13
  res$funvars[, "Heritability"] <- c(0.71741, 0.62372, 0.70984)
  out <- capture.output(print(summary(res)))

  ## the note names the offending function
  expect_true(any(grepl("Note: for Heritability", out)))

  ## the note must say *why*, not merely flag the gap
  expect_true(any(grepl("non-linear", out)))
  expect_true(any(grepl("boundary", out)))
})
