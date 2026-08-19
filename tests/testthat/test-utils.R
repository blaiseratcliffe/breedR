### Test the auxiliar functions in utils.R ###


context("Auxiliar functions")


test_that("lmat2df() works as expected", {
  
  rnms <- letters[1:3]
  cnms <- c("x", "y")
  label <- "e"
  
  tm <- matrix(1:9, 3, 3, dimnames = rep(list(rnms), 2))
  
  tlm <- structure(rep(list(tm), 2), names = cnms)
  
  tdf <- lmat2df(tlm, label)

  ## labels
  exp_labs <- c("e.a", "e.a_e.b", "e.a_e.c", "e.b", "e.b_e.c", "e.c")
  expect_identical(rownames(tdf), exp_labs)
  
  ## columns
  expect_identical(colnames(tdf), cnms)
  
  ## values
  ## only lower-triangular values
  exp_values <- tm[lower.tri(tm, diag = TRUE)]
  expect_identical(tdf$x, exp_values)
  expect_identical(tdf$y, exp_values)
})


#### Working directories ####

test_that("breedR_workdir() gives a fresh directory under tempdir()", {

  d1 <- breedR_workdir()
  d2 <- breedR_workdir()
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)

  expect_true(dir.exists(d1))
  expect_true(dir.exists(d2))

  ## Distinct: the whole point is that a second run cannot land on the first
  ## one's files.
  expect_false(identical(d1, d2))

  expect_true(startsWith(normalizePath(d1, winslash = "/"),
                         normalizePath(tempdir(), winslash = "/")))

  ## the prefix reaches the name, so a stray directory can be traced back
  expect_match(basename(breedR_workdir('breedR_gibbs_')), '^breedR_gibbs_')
})


test_that("clean_workdir() accepts every shape the package produces", {

  ## a bare path
  d <- breedR_workdir()
  expect_true(clean_workdir(d))
  expect_false(dir.exists(d))

  ## a fitted model
  d <- breedR_workdir()
  expect_true(clean_workdir(list(reml = list(dir = d))))
  expect_false(dir.exists(d))

  ## a gibbsf90()/postgsf90()/renumf90() result
  d <- breedR_workdir()
  expect_true(clean_workdir(list(dir = d)))
  expect_false(dir.exists(d))
})


test_that("clean_workdir() is a no-op when there is nothing to remove", {

  ## already gone
  d <- breedR_workdir()
  unlink(d, recursive = TRUE)
  expect_false(clean_workdir(d))

  ## objects that never recorded a directory -- e.g. anything rebuilt from a
  ## saved file by an older version of the package
  expect_false(clean_workdir(list()))
  expect_false(clean_workdir(list(reml = list())))
})


test_that("clean_workdir() refuses to leave tempdir()", {

  ## The guard that matters. renumf90() and gibbsf90_from_renum() take a
  ## user-chosen `dir`, so a stored path is not automatically ours to delete.
  outside <- file.path(getwd(), 'clean_workdir_outside')
  dir.create(outside, showWarnings = FALSE)
  on.exit(unlink(outside, recursive = TRUE), add = TRUE)
  writeLines('keep me', file.path(outside, 'data.txt'))

  expect_error(clean_workdir(outside), 'Refusing to remove')
  expect_true(file.exists(file.path(outside, 'data.txt')))

  ## and tempdir() itself is not a working directory either -- removing it
  ## would take every other fit in the session with it
  expect_error(clean_workdir(tempdir()), 'Refusing to remove')
  expect_true(dir.exists(tempdir()))
})
