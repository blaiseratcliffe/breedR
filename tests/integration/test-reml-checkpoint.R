context("REML progress and resume")

## These require the PROGSF90 binaries.

dat <- breedR::globulus
ped <- list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self')

fit_globulus <- function(...)
  suppressMessages(
    remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped, data = dat, ...)
  )

## A run whose wall-clock and CPU stamps differ from every other run's.
volatile_line <- function(x)
  grepl("^ \\*|in +[0-9.]+ +s,|CPU TIME|TIME FOR |elapsed time", x)

stable <- function(x) x[!volatile_line(x)]

wd <- file.path(tempdir(), "reml_checkpoint")
unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE)
## remlf90() returns the path normalised with forward slashes
wd <- normalizePath(wd, winslash = "/")

## Load the same copy of breedR in a spawned R process as the one under test:
## during development that is the source tree, in CI the installed package.
loader <- local({
  src <- tryCatch(if (requireNamespace("pkgload", quietly = TRUE))
                    pkgload::pkg_path() else NULL,
                  error = function(e) NULL)
  if (!is.null(src))
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))", shQuote(src))
  else "suppressMessages(library(breedR))"
})


test_that("progress_file streams the backend output to disk", {

  f <- file.path(wd, "stream.log")
  res <- fit_globulus(progress_file = f)

  expect_true(file.exists(f))
  expect_identical(res$reml$progress_file, f)

  logged <- readLines(f, warn = FALSE)
  expect_true(any(grepl("In round", logged)))
  expect_true(any(grepl("Final Estimates", logged)))
  expect_identical(res$reml$output, logged)
})


test_that("streaming does not change what the parser sees", {

  ## The interesting assertion is not that reml.out equals the file it was
  ## written from -- that is true by construction -- but that a streamed run
  ## yields the same output as an ordinary one. This is what would catch
  ## stderr being folded into the captured stream.
  plain <- fit_globulus()
  streamed <- fit_globulus(progress_file = file.path(wd, "cmp.log"))
  plain2 <- fit_globulus()

  ## Sanity: the filter removes exactly the lines that are unstable anyway.
  expect_identical(stable(plain$reml$output), stable(plain2$reml$output))

  expect_identical(stable(plain$reml$output), stable(streamed$reml$output))
  expect_equal(plain$reml$rounds, streamed$reml$rounds)
  expect_equal(plain$var, streamed$var)
})


test_that("the log is readable while the fit is still running", {

  ## The regression test for the Windows behaviour of system2(stdout = <file>),
  ## which holds an exclusive lock: readLines(), file.copy() and `type` all
  ## fail with permission denied until the process exits. Tailing a multi-day
  ## fit is the point of the argument, so this has to keep working.
  f <- file.path(wd, "live.log")
  script <- file.path(wd, "live_fit.R")

  writeLines(c(
    loader,
    "dat <- breedR::globulus",
    "invisible(suppressMessages(remlf90(",
    "  fixed = phe_X ~ gg, random = ~ bl,",
    "  genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),",
    "  spatial = list(model = 'AR', coord = dat[, c('x','y')], rho = c(.85,.8)),",
    "  data = dat, method = 'em',",
    sprintf("  progress_file = %s)))", shQuote(f))
  ), script)

  unlink(f)
  spawn_log <- file.path(wd, "live_fit.out")
  system2(file.path(R.home("bin"), "Rscript"), args = shQuote(script),
          wait = FALSE, stdout = spawn_log, stderr = spawn_log)

  rounds_seen <- 0L
  finished <- FALSE
  deadline <- Sys.time() + 120

  while (Sys.time() < deadline) {
    Sys.sleep(0.1)
    if (!file.exists(f)) next
    ## The assertion: this read must not fail while another process writes.
    seen <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
    expect_false(is.null(seen))
    if (is.null(seen)) break
    finished <- any(grepl("Final Estimates", seen))
    if (!finished) rounds_seen <- max(rounds_seen, sum(grepl("In round", seen)))
    if (finished) break
  }

  ## Surface the spawned process's own diagnostics if it never got going.
  if (!finished && file.exists(spawn_log))
    message("spawned fit said: ",
            paste(utils::tail(readLines(spawn_log, warn = FALSE), 10),
                  collapse = " | "))

  expect_true(finished)
  expect_gt(rounds_seen, 0L)
})


test_that("a relative progress_file lands in the calling directory", {

  ## remlf90() moves into tempdir() before running the backend, so a path
  ## resolved too late would be written there and lost with the session.
  d <- file.path(wd, "relative")
  dir.create(d, showWarnings = FALSE)
  owd <- setwd(d)
  on.exit(setwd(owd))

  res <- fit_globulus(progress_file = "rel.log")

  expect_true(file.exists(file.path(d, "rel.log")))
  expect_false(file.exists(file.path(tempdir(), "rel.log")))
  expect_identical(basename(res$reml$progress_file), "rel.log")
})


test_that("reml_checkpoint() recovers the reported variance components", {

  f <- file.path(wd, "ckpt.log")
  res <- fit_globulus(progress_file = f)

  vc <- reml_checkpoint(f, model = res)

  expect_identical(names(vc), c("bl", "genetic", "residuals"))
  expect_equal(attr(vc, "round"), res$reml$rounds)
  expect_equal(unname(vapply(vc, function(m) m[1, 1], 1)),
               unname(res$var[, "Estimated variances"]),
               tolerance = 1e-6)

  ## and without a model, positionally
  vp <- reml_checkpoint(f)
  expect_identical(names(vp), c("G1", "G2", "residuals"))
})


test_that("cont resumes from the previous log", {

  f <- file.path(wd, "resume.log")
  cold <- fit_globulus(progress_file = f)

  expect_message(
    warm <- remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped,
                    data = dat, progress_file = f, cont = TRUE),
    "Resuming from round"
  )

  expect_equal(warm$reml$resumed_from, cold$reml$rounds)
  expect_lt(warm$reml$rounds, cold$reml$rounds)
  expect_equal(warm$var[, "Estimated variances"],
               cold$var[, "Estimated variances"],
               tolerance = 1e-4)

  ## the previous log is preserved rather than clobbered
  expect_true(any(grepl("^resume\\.log\\.[0-9]{8}-[0-9]{6}$",
                        list.files(wd))))

  ## and the recovered residual reaches the parameter file
  pars <- readLines(file.path(tempdir(), "parameters"), warn = FALSE)
  expect_true(any(grepl("RANDOM_RESIDUAL VALUES", pars)))
})


test_that("cont works after a non-converged run", {

  ## The best case for resuming: hitting the iteration cap returns all-NA
  ## variance components, while the log holds every round intact.
  f <- file.path(wd, "em.log")
  cold <- fit_globulus(progress_file = f, method = 'em',
                       progsf90.options = 'maxrounds 12')

  expect_message(
    warm <- remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped,
                    data = dat, method = 'em',
                    progress_file = f, cont = TRUE),
    "Resuming from round"
  )
  expect_true(warm$reml$resumed_from > 0L)
})


test_that("the progress and resume arguments are guarded", {

  f <- file.path(wd, "guard.log")
  fit_globulus(progress_file = f)

  ## cont together with an explicit (complete) var.ini
  expect_error(
    suppressMessages(
      remlf90(fixed = phe_X ~ gg, random = ~ bl,
              genetic = c(ped, list(var.ini = 1)), data = dat,
              var.ini = list(bl = 1, residuals = 1),
              progress_file = f, cont = TRUE)),
    "cont = TRUE"
  )

  expect_error(fit_globulus(progress_file = file.path(wd, "nope.log"),
                            cont = TRUE),
               "does not exist")

  expect_error(fit_globulus(cont = TRUE), "requires 'progress_file'")

  ## an AR rho grid is N fits, not one
  expect_error(
    suppressMessages(
      remlf90(fixed = phe_X ~ gg, data = dat,
              spatial = list(model = 'AR', coord = dat[, c('x', 'y')],
                             rho = rbind(c(.8, .8), c(.9, .9))),
              progress_file = file.path(wd, "grid.log"))),
    "rho grid"
  )
})


test_that("a non-zero exit from the streamed backend is surfaced", {

  ## When streaming, the exit status comes back from close(<pipe>) rather than
  ## as an attribute of the captured output, so the plumbing differs from the
  ## plain path and needs its own cover.
  ##
  ## It cannot be exercised through remlf90(): BLUPF90+ exits 0 even on fatal
  ## input errors ("There is no such data file", "NUMBER_OF_EFFECTS not
  ## found"), which is why breedR reports those as parse failures instead.
  ## Only an outright crash returns non-zero. So drive the mechanism directly.
  rs <- file.path(R.home("bin"), "Rscript")
  d <- file.path(wd, "status"); dir.create(d, showWarnings = FALSE)
  owd <- setwd(d); on.exit(setwd(owd))

  writeLines('cat("some output\\n"); quit(status = 7)', "boom.R")
  writeLines("parameters", "pf90_stdin")

  con <- pipe(paste0(shQuote(rs), ' boom.R < pf90_stdin 2> pf90_stderr'), 'r')
  out <- readLines(con, warn = FALSE)
  status <- close(con)

  expect_equal(as.integer(status), 7L)

  attr(out, 'status') <- as.integer(status)
  expect_error(stop_progsf90_failure(out, attr(out, 'status')),
               "exit code 7")
  expect_error(stop_progsf90_failure(out, attr(out, 'status')),
               "some output")
})
