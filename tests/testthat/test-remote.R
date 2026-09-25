### Remote execution plumbing ###

context("Remote computing")

## The remote path needs a configured Linux server, so nothing here fits a
## model. What can be checked without one is the property the fix relies on.

test_that("results are never retrieved into a shared directory", {

  ## retrieve_remote() used to default `dest` to tempdir(), and breedR.qget()
  ## took that default: every job of a session was untarred into the same
  ## place, so a second retrieval overwrote the first one's LOG and solutions
  ## and both objects then advertised the one path. Every caller now passes a
  ## directory of its own, and `dest` has no default so that a third caller
  ## cannot silently go back to sharing one.
  expect_true(identical(formals(retrieve_remote)$dest, quote(expr = )))
  expect_true(identical(formals(breedR.remote)$dest, quote(expr = )))

  ## and the one caller that has no per-fit directory to pass makes itself one
  expect_match(paste(deparse(body(breedR.qget)), collapse = ' '),
               'retrieve_remote\\(rdir, dest = breedR_workdir')
})


test_that("breedR.qget() surfaces the backend's own diagnostic when a submitted run never started (issue #14)", {

  ## A remote fit that dies during setup (#14) is reported "Finished" by
  ## qstat -- the backend exited 0 -- but never wrote a solutions file.
  ## breedR.qget() used to call parse_results() unguarded and fail on
  ## read.table()'s "cannot open the connection" instead of the backend's own
  ## diagnostic, which is sitting in the retrieved LOG.
  fake_dir <- file.path(tempdir(), "breedR_qget_test_14")
  dir.create(fake_dir, showWarnings = FALSE)
  withr::defer(unlink(fake_dir, recursive = TRUE))
  writeLines(c("some setup echo",
              "There is no such data file: no_such_file"),
            file.path(fake_dir, "LOG"))
  # deliberately no 'solutions' file

  data <- data.frame(x = seq(0, 2, length.out = 100),
                     g = factor(rep(1:50, 2)),
                     y = rep(c(9, 11), 50))
  mc <- call('remlf90', fixed = quote(y ~ x), random = quote(~ g), data = quote(data))
  mf <- build.mf(mc)
  effects <- build.effects(mf, NULL, NULL, NULL, list(g = 3.4))

  local_mocked_bindings(
    breedR.qstat = function(id) list(list(id = "1", status = "Finished", pid = "123")),
    retrieve_remote = function(rdir, dest) fake_dir,
    breedR.qdel = function(id) invisible(NULL),
    .package = 'breedR'
  )

  fake_id <- structure(list(id = "1", effects = effects, mf = mf,
                            method = 'ai', mcout = quote(remlf90())),
                       class = c('breedR', 'remlf90'))

  expect_error(breedR.qget(fake_id), "no such data file", fixed = TRUE)
})


test_that("breedR.qget() applies the same non-convergence check to a retrieved log (issue #40)", {

  ## A remote fit stopped by its iteration cap (#40) is as unmarked in the
  ## retrieved LOG as in a local one. breedR.qget() has only that log, not
  ## the options of the fit, and must still warn and give no estimates.
  fake_dir <- file.path(tempdir(), "breedR_qget_test_40")
  dir.create(fake_dir, showWarnings = FALSE)
  withr::defer(unlink(fake_dir, recursive = TRUE))
  file.copy(file.path(testdata, "issue40_capped_ai.log"),
            file.path(fake_dir, "LOG"))
  file.copy(file.path(testdata, "issue40_capped_ai.sol"),
            file.path(fake_dir, "solutions"))

  data <- issue40_data()
  mc <- call('remlf90', fixed = quote(y ~ x), random = quote(~ g), data = quote(data))
  mf <- build.mf(mc)
  effects <- build.effects(mf, NULL, NULL, NULL, list(g = 3.4))

  local_mocked_bindings(
    breedR.qstat = function(id) list(list(id = "1", status = "Finished", pid = "123")),
    retrieve_remote = function(rdir, dest) fake_dir,
    breedR.qdel = function(id) invisible(NULL),
    .package = 'breedR'
  )

  fake_id <- structure(list(id = "1", effects = effects, mf = mf,
                            method = 'ai', mcout = quote(remlf90())),
                       class = c('breedR', 'remlf90'))

  expect_warning(res <- suppressMessages(breedR.qget(fake_id)),
                 "did not converge")
  expect_true(all(is.na(res$var)))
})
