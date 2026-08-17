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

  ## remlf90() moves into its working directory before running the backend, so
  ## a path resolved too late would be written there and lost with the session.
  d <- file.path(wd, "relative")
  dir.create(d, showWarnings = FALSE)
  owd <- setwd(d)
  on.exit(setwd(owd))

  res <- fit_globulus(progress_file = "rel.log")

  expect_true(file.exists(file.path(d, "rel.log")))
  ## against the fit's own directory, not tempdir(): that is where a late
  ## resolution would now land it, so checking tempdir() would pass whether or
  ## not the bug came back.
  expect_false(file.exists(file.path(res$reml$dir, "rel.log")))
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

  ## and the recovered residual reaches the parameter file of the resumed fit
  pars <- readLines(file.path(warm$reml$dir, "parameters"), warn = FALSE)
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

  ## When streaming, the exit status comes back from the process object rather
  ## than as an attribute of the captured output, so the plumbing differs from
  ## the plain path and needs its own cover.
  ##
  ## It cannot be exercised through remlf90(): BLUPF90+ exits 0 even on fatal
  ## input errors ("There is no such data file", "NUMBER_OF_EFFECTS not
  ## found"), which is why breedR reports those as parse failures instead.
  ## Only an outright crash returns non-zero. So drive the mechanism directly,
  ## in the same shape remlf90() uses it.
  rs <- file.path(R.home("bin"), "Rscript")
  d <- file.path(wd, "status"); dir.create(d, showWarnings = FALSE)
  owd <- setwd(d); on.exit(setwd(owd))

  writeLines('cat("some output\\n"); quit(status = 7)', "boom.R")
  writeLines("parameters", "pf90_stdin")

  px <- processx::process$new(rs, "boom.R", stdin = "pf90_stdin",
                              stdout = "|", stderr = "pf90_stderr")
  out <- character(0)
  repeat {
    px$poll_io(1000)
    l <- px$read_output_lines()
    if (length(l)) out <- c(out, l) else if (!px$is_alive()) break
  }
  status <- px$get_exit_status()

  expect_equal(as.integer(status), 7L)

  attr(out, 'status') <- as.integer(status)
  expect_error(stop_progsf90_failure(out, attr(out, 'status')),
               "exit code 7")
  expect_error(stop_progsf90_failure(out, attr(out, 'status')),
               "some output")
})


test_that("a failure to open the log does not leave a backend running", {

  ## The log is opened before the backend is started, so the one fallible step
  ## in that sequence cannot orphan a child. An orphan holds files in
  ## tempdir(), and the symptom is that the *next* fit in the same session
  ## fails with "The process cannot access the file because it is being used
  ## by another process" -- so that is what this asserts.
  d <- file.path(wd, "not_a_file.log")
  dir.create(d, showWarnings = FALSE)

  ## file() also warns; the error is what matters here
  expect_error(suppressWarnings(fit_globulus(progress_file = d)))

  ## the session must still be usable
  expect_error(after <- fit_globulus(), NA)
  expect_true(after$reml$rounds > 0L)
})


test_that("a failed resume leaves the previous log where it was", {

  ## The log is the artefact the whole feature exists to protect, so nothing
  ## that can still reject the run may run after it has been moved aside.
  f <- file.path(wd, "survive.log")
  fit_globulus(progress_file = f)
  expect_true(file.exists(f))

  ## resume with a model whose random groups differ from the log's
  expect_error(
    suppressMessages(
      remlf90(fixed = phe_X ~ gg, genetic = ped, data = dat,
              progress_file = f, cont = TRUE)),
    "Cannot resume"
  )

  ## the error named this file; it must still be there
  expect_true(file.exists(f))
  expect_false(any(grepl("^survive[.]log[.][0-9]{8}-[0-9]{6}$", list.files(wd))))

  ## and the retry the user will actually type must work
  expect_message(ok <- remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped,
                               data = dat, progress_file = f, cont = TRUE),
                 "Resuming from round")
  expect_true(ok$reml$resumed_from > 0L)
})


test_that("a resume that cannot open its log puts the previous one back", {

  ## The rename is the last step before the log is opened, and that open can
  ## still fail; without a restore the user would be left with a timestamped
  ## backup and no file under the name their retry uses.
  f <- file.path(wd, "restore.log")
  fit_globulus(progress_file = f)
  before <- readLines(f, warn = FALSE)

  ## make the open fail: hold the path open elsewhere is unreliable across
  ## platforms, so drop a directory in the way after the log has been written
  lock <- file.path(wd, "restore_lock")
  dir.create(lock, showWarnings = FALSE)

  expect_error(suppressWarnings(
    fit_globulus(progress_file = lock, cont = TRUE)))

  ## the real log is untouched, and its own resume still works
  expect_true(file.exists(f))
  expect_identical(readLines(f, warn = FALSE), before)
})


test_that("cont preserves the previous run's error stream too", {

  ## .err normally holds the diagnostic that prompted the resume, so keeping
  ## the log while deleting it would preserve the wrong half.
  f <- file.path(wd, "err.log")
  fit_globulus(progress_file = f)
  writeLines("forrtl: severe: something went wrong", paste0(f, ".err"))

  expect_message(remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped,
                         data = dat, progress_file = f, cont = TRUE),
                 "Resuming from round")

  bak <- grep("^err[.]log[.][0-9]{8}-[0-9]{6}$", list.files(wd), value = TRUE)
  expect_length(bak, 1L)
  expect_true(file.exists(file.path(wd, paste0(bak, ".err"))))
  expect_match(readLines(file.path(wd, paste0(bak, ".err")))[1], "forrtl")
})


test_that("a resume with no error stream warns about nothing", {

  ## file.rename() on a missing source returns FALSE *and* warns. Most resumes
  ## have no .err -- the documented flagship case, a fit that exhausted its
  ## iterations, emits no stderr at all -- so an unguarded rename would warn
  ## on nearly every resume.
  f <- file.path(wd, "noerr.log")
  fit_globulus(progress_file = f)
  unlink(paste0(f, ".err"))

  expect_warning(
    suppressMessages(
      remlf90(fixed = phe_X ~ gg, random = ~ bl, genetic = ped,
              data = dat, progress_file = f, cont = TRUE)),
    NA)
})


test_that("killing the backend frees the session it was poisoning", {

  ## Abandoning a running backend is not benign: once its output buffer fills
  ## with nobody reading, it blocks forever holding that fit's working files,
  ## and a later run against them fails with "Permission denied" on
  ## parameters. That is why remlf90() registers a kill on exit.
  ##
  ## Scope: this drives the mechanism rather than remlf90(), because remlf90()
  ## exposes no handle on the process and there is no portable way to deliver
  ## an interrupt to our own R session from testthat. What it does reproduce
  ## faithfully is the poisoning and its cure, in one session, against the real
  ## backend -- which is the part that used to be broken.
  ##
  ## Run in the priming fit's own directory: that is where its parameter file
  ## is, and the backend below is started on it by name.
  prime <- fit_globulus(progress_file = file.path(wd, "prime.log"))
  td <- prime$reml$dir
  owd <- setwd(td); on.exit(setwd(owd), add = TRUE)
  writeLines("parameters", "pf90_stdin")

  bin <- file.path(breedR.getOption("breedR.bin"),
                   progsf90_files(breedR.os.type()))
  px <- processx::process$new(bin, stdin = "pf90_stdin",
                              stdout = "|", stderr = "pf90_stderr")

  ## read a little, then walk away exactly as an interrupt would
  seen <- 0L
  deadline <- Sys.time() + 60
  repeat {
    px$poll_io(1000)
    l <- px$read_output_lines()
    seen <- seen + sum(grepl("In round", l))
    if (seen >= 3L || !px$is_alive() || Sys.time() > deadline) break
  }
  expect_true(px$is_alive())

  t0 <- Sys.time()
  px$kill()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  expect_lt(elapsed, 5)            # must not wait for the fit to finish
  expect_false(px$is_alive())

  ## the files it held are released ...
  expect_error({ h <- file(file.path(td, "parameters"), "w"); close(h) }, NA)

  ## ... and the session can still fit, which is the symptom users hit
  expect_error(after <- fit_globulus(), NA)
  expect_true(after$reml$rounds > 0L)
})
