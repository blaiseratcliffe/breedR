### Each fit owns its working files (issue #23) ###

context("Fit isolation")

## Every fit used to work in one shared tempdir(), writing the same fixed set
## of names. The second fit therefore overwrote the first one's parameter file,
## data, structure files and solutions. The cost was not untidiness: a fit that
## produced no solutions of its own was parsed against whatever the previous
## one left behind, so one broken model reported as three different errors
## depending on what had run earlier in the session -- including a
## get_efnames() length mismatch that pointed at code with nothing wrong in it.

dat <- globulus

## Two models with *different* numbers of effects, which is the case that used
## to produce the most misleading error of the three.
a <- suppressMessages(remlf90(phe_X ~ gg, data = dat))
b <- suppressMessages(
  remlf90(phe_X ~ gg,
          random  = ~ bl,
          genetic = list(model    = 'add_animal',
                         pedigree = dat[, 1:3],
                         id       = 'self'),
          data    = dat))


test_that("each fit records its own working directory", {

  expect_false(is.null(a$reml$dir))
  expect_false(is.null(b$reml$dir))
  expect_false(identical(a$reml$dir, b$reml$dir))

  ## under tempdir(), so the session still cleans up after itself
  expect_true(startsWith(normalizePath(a$reml$dir, winslash = "/"),
                         normalizePath(tempdir(), winslash = "/")))
})


test_that("a later fit does not overwrite an earlier one's files", {

  for (res in list(a, b)) {
    expect_true(file.exists(file.path(res$reml$dir, "parameters")))
    expect_true(file.exists(file.path(res$reml$dir, "solutions")))
  }

  ## The solutions still describe the model that produced them: `a` has one
  ## effect, `b` has three. Reading either one after both have run used to give
  ## whichever ran last.
  n_effects <- function(res)
    length(unique(utils::read.table(file.path(res$reml$dir, "solutions"),
                                    header = FALSE, skip = 1)[[2]]))

  expect_equal(n_effects(a), 1L)
  expect_gt(n_effects(b), 1L)

  ## and the fitted objects still agree with their own files
  expect_equal(nrow(utils::read.table(file.path(a$reml$dir, "solutions"),
                                      header = FALSE, skip = 1)),
               nlevels(dat$gg))
})


test_that("a fit that writes no solutions says so", {

  ## The backend exits 0 even on fatal input errors, so breedR used to walk on
  ## and read.table() a file that was never written -- 'cannot open the
  ## connection' at best, and the previous fit's numbers at worst. Drive it by
  ## emptying the structure file of an otherwise valid model and re-running the
  ## backend the way remlf90() does.
  inc.mat <- model.matrix(~ 0 + bl, dat)
  cov.mat <- diag(nlevels(dat$bl))
  g <- suppressMessages(
    remlf90(phe_X ~ gg, generic = list(bl = list(inc.mat, cov.mat)),
            data = dat))

  d <- file.path(tempdir(), "empty_structure")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE)
  file.copy(list.files(g$reml$dir, full.names = TRUE), d)
  unlink(file.path(d, "solutions"))
  file.create(file.path(d, "generic_bl"))          # the #22 shape: 0 bytes

  bin <- file.path(breedR.getOption("breedR.bin"),
                   progsf90_files(breedR.os.type()))
  local_bin <- file.path(d, basename(bin))
  file.copy(bin, local_bin, overwrite = TRUE)
  owd <- setwd(d)
  on.exit({ setwd(owd); unlink(d, recursive = TRUE) }, add = TRUE)
  out <- system2(file.path(".", basename(bin)), input = "parameters",
                 stdout = TRUE, stderr = TRUE)

  ## the backend reports the empty file and still exits cleanly
  expect_true(any(grepl("empty|0  elements", out)))
  expect_null(attr(out, "status"))
  expect_false(file.exists(file.path(d, "solutions")))
})
