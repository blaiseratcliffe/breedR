context("REML checkpointing")

log_lines <- function(f) readLines(file.path(testdata, f), warn = FALSE)

x1 <- log_lines("airemlf90_log_1.txt")   # 2 traits, no random effects
x3 <- log_lines("airemlf90_log_3.txt")   # 2 traits, 3 random groups
x4 <- log_lines("airemlf90_log_4.txt")   # 10 traits, wrapped rows
y1 <- log_lines("remlf90_log_1.txt")     # legacy REMLF90, EM dialect


## -- round anchors --------------------------------------------------------

test_that("reml_round_index() finds every round header", {

  expect_error(idx <- reml_round_index(x3), NA)
  expect_equal(length(idx), 18L)
  expect_identical(names(idx), as.character(1:18))
  expect_true(all(grepl("In round", x3[idx])))

  expect_equal(length(reml_round_index(x1)), 8L)
  expect_equal(length(reml_round_index(x4)), 6L)

  ## the legacy log prints both " In round N" and " round N"; anchoring on
  ## the former must still give one segment per round
  expect_equal(length(reml_round_index(y1)), 542L)

  expect_identical(reml_round_index(character(0)), structure(integer(0),
                                                             names = character(0)))
})


## -- the corruption guard -------------------------------------------------

test_that("a block truncated inside an exponent is rejected, not parsed", {

  ## This is the failure mode a crash actually produces: the log ends in the
  ## middle of a line. parse.txtmat() does not error on it -- ' 0.57309E+09'
  ## cut to ' 0.573' parses cleanly and is wrong by nine orders of magnitude.
  mk <- grep("^[[:space:]]*(new [Gg]|G)[[:space:]]*$", x3)
  rr <- grep("In round", x3)
  last <- tail(mk[mk > tail(rr, 1)], 1)

  truncated <- c(x3[1:(last + 1)], substr(x3[last + 2], 1, 20))

  ## the raw parser is happy to return garbage ...
  raw <- parse.txtmat(extract_block(last + 1L, truncated))
  expect_true(is.matrix(raw))
  expect_equal(raw[2, 2], 0.573)          # should have been 5.7309e8

  ## ... the strict wrapper refuses it
  expect_null(parse_block_strict(last + 1L, truncated))

  ## and the round is skipped rather than trusted
  ck <- last_reml_checkpoint(truncated, n_groups = 3L, spd = FALSE)
  expect_true(attr(ck, "round") < 18L)
})


test_that("parse_block_strict() refuses any block running to end of file", {

  ## A finished block is always followed by more output. Refusing blocks that
  ## reach the last line is what makes the parser safe on a log still being
  ## written, where the tail is incomplete by construction.
  mk <- grep("^[[:space:]]*new G[[:space:]]*$", x3)

  expect_null(parse_block_strict(mk[1] + 1L, x3[1:(mk[1] + 2)]))
  expect_true(is.matrix(parse_block_strict(mk[1] + 1L, x3[1:(mk[1] + 3)])))
})


test_that("both clean truncation modes are caught", {

  mk <- grep("^[[:space:]]*new G[[:space:]]*$", x3)[1]

  ## one row of a 2x2, followed by more text: parse.txtmat() errors
  one_row <- c(x3[1:(mk + 1)], " new G", "  1.0  2.0", "  2.0  5.0", "end")
  expect_error(parse.txtmat(extract_block(mk + 1L, one_row)),
               "square matrix dimensions")
  expect_null(parse_block_strict(mk + 1L, one_row))

  ## heading with no rows at all: extract_block() errors
  no_rows <- c(x3[1:mk], "some text", "more text")
  expect_error(extract_block(mk + 1L, no_rows), "is_numericlog")
  expect_null(parse_block_strict(mk + 1L, no_rows))
})


## -- the round-trip invariant ---------------------------------------------

test_that("the last round equals the final estimates", {

  ## This is what makes resuming possible at all: what the optimiser prints
  ## each round is exactly what it reports at the end.
  ck <- last_reml_checkpoint(x3, n_groups = 3L, spd = FALSE)
  expect_equal(attr(ck, "round"), 18L)

  final <- lapply(grep("Genetic variance|Residual variance", x3) + 1L,
                  function(i) parse.txtmat(extract_block(i, x3)))

  expect_equal(length(ck), length(final))
  for (i in seq_along(final))
    expect_equal(unname(ck[[i]]), unname(final[[i]]))
})


test_that("the residual is placed last, though it is printed first", {

  ## Within a round the order is `new R` then the `new G`s; in the final
  ## estimates the residual comes last. The two are not interchangeable.
  ck <- last_reml_checkpoint(x3, n_groups = 3L, spd = FALSE,
                             names = c("a", "b", "c", "residuals"))
  expect_identical(names(ck), c("a", "b", "c", "residuals"))

  resid_final <- parse.txtmat(extract_block(grep("Residual variance", x3) + 1L, x3))
  expect_equal(unname(ck$residuals), unname(resid_final))
})


## -- dialects -------------------------------------------------------------

test_that("all three round-block dialects are read", {

  ## BLUPF90+ under AI-REML: ' new R' / ' new G'
  expect_equal(attr(last_reml_checkpoint(x3, 3L, spd = FALSE), "round"), 18L)

  ## Legacy REMLF90 under EM: ' new r' / bare ' G'
  expect_equal(attr(last_reml_checkpoint(y1, 1L, spd = FALSE), "round"), 541L)

  ## Current BLUPF90+ under EM-REML mixes them: ' new r' with ' new G'.
  ## No fixture shows this, and a parser matching only ' new R'/' new G'
  ## silently recovers nothing from a method = 'em' run.
  mixed <- c(" In round            1  convergence=  1.0",
             " delta convergence=  1.0",
             " new r", "  9.0",
             " new G", "  1.0",
             " -2logL =    1.0       : AIC =    1.0",
             " In round            2  convergence=  0.1",
             " delta convergence=  0.1",
             " new r", "  8.0",
             " new G", "  2.0",
             " solutions stored")

  ck <- last_reml_checkpoint(mixed, n_groups = 1L)
  expect_equal(attr(ck, "round"), 2L)
  expect_equal(unname(ck[[1]]), matrix(2))
  expect_equal(unname(ck$residuals), matrix(8))
})


test_that("a bare G heading is not confused with a diagnostic message", {

  ## ' G not positive definite: fixed (setup_g)' appears in real logs and
  ## must not be taken for the legacy bare-G block heading.
  noise <- grep("not positive definite", x3, value = TRUE)
  expect_true(length(noise) > 0)
  expect_false(any(grepl("^[[:space:]]*(new [Gg]|G)[[:space:]]*$", noise)))
})


## -- positive definiteness ------------------------------------------------

test_that("non positive definite rounds are skipped", {

  ## Rounds 13-18 of this fixture are all non-SPD, so the walk back from a
  ## converged fit lands six rounds earlier.
  expect_equal(attr(last_reml_checkpoint(x3, 3L, spd = FALSE), "round"), 18L)
  expect_equal(attr(last_reml_checkpoint(x3, 3L, spd = TRUE), "round"), 12L)
})


## -- shapes ---------------------------------------------------------------

test_that("a model with no random effects yields the residual alone", {

  ck <- last_reml_checkpoint(x1, n_groups = 0L)
  expect_equal(attr(ck, "round"), 8L)
  expect_identical(names(ck), "residuals")
  expect_equal(unname(ck$residuals),
               matrix(c(0.94826, 2.8013, 2.8013, 12.120), 2))
})


test_that("scalar blocks give 1x1 matrices", {

  ## No fixture is univariate with a single random effect, but this is the
  ## commonest real model.
  uni <- c(" In round            1  convergence=  1.0",
           " delta convergence=  1.0",
           " new R", "  10.5",
           " new G", "  2.66",
           " -2logL =    1.0       : AIC =    1.0",
           " In round            2  convergence=  0.1",
           " delta convergence=  0.1",
           " new R", "  9.5",
           " new G", "  3.10",
           " solutions stored")

  ck <- last_reml_checkpoint(uni, n_groups = 1L, dims = c(1L, 1L))
  expect_equal(attr(ck, "round"), 2L)
  expect_equal(dim(ck[[1]]), c(1L, 1L))
  expect_equal(unname(ck[[1]]), matrix(3.10))
  expect_equal(unname(ck$residuals), matrix(9.5))
})


test_that("wrapped matrix rows are reassembled", {

  ## 10 traits: each row of the residual matrix spans two physical lines.
  ck <- last_reml_checkpoint(x4, n_groups = 0L, spd = FALSE)
  expect_equal(attr(ck, "round"), 6L)
  expect_equal(dim(ck$residuals), c(10L, 10L))
})


## -- refusals -------------------------------------------------------------

test_that("an unusable log yields an empty result and a reason", {

  ## NULL cannot carry attributes, so failure is a zero-length list.
  none <- last_reml_checkpoint(character(0), 0L)
  expect_equal(length(none), 0L)
  expect_identical(attr(none, "reason"), "no completed round")

  ## round 1 of this fixture starts at line 103
  early <- last_reml_checkpoint(head(x3, 100), 3L)
  expect_equal(length(early), 0L)

  ## a log written by a model with a different number of groups
  wrong <- last_reml_checkpoint(x3, n_groups = 2L, spd = FALSE)
  expect_equal(length(wrong), 0L)
  expect_identical(attr(wrong, "reason"), "dimension")
})


test_that("the group count is not inferred from a truncated round", {

  ## Without a model, reml_checkpoint() has to work out how many groups to
  ## expect. Reading that off the last round is wrong: a round cut short by a
  ## crash shows too few, and the result is a checkpoint quietly missing a
  ## component rather than a fall-back to the previous round.
  full <- c(" In round            1  convergence=  1.0",
            " new R", "  9.0",
            " new G", "  1.0",
            " new G", "  2.0",
            " -2logL =    1.0       : AIC =    1.0",
            " In round            2  convergence=  0.1",
            " new R", "  8.0",
            " new G", "  3.0",
            " new G", "  4.0",
            " solutions stored")

  ck <- last_reml_checkpoint(full, n_groups = 2L)
  expect_equal(attr(ck, "round"), 2L)

  ## now cut round 2 after its first group
  cut <- c(head(full, 13), "  9.9")
  f <- tempfile(fileext = ".log")
  writeLines(cut, f)
  on.exit(unlink(f))

  got <- reml_checkpoint(f)
  expect_equal(length(got), 3L)             # two groups plus the residual
  expect_equal(attr(got, "round"), 1L)      # not round 2, which lost a group
})


test_that("a log from a different model is refused by name", {

  ## one random group, but the log records three
  eff <- list(bl = effect_group(list(diagonal(factor(rep(1:3, 4)))),
                                cov.ini = matrix(1), ntraits = 1L))
  expect_error(
    reml_checkpoint_var.ini("somewhere.log", eff, ntraits = 1L, lines = x3),
    "different model")
})


## -- argument guards ------------------------------------------------------

test_that("check_progress_args() rejects bad combinations", {

  expect_error(check_progress_args(NULL, TRUE, "/bin"),
               "requires 'progress_file'")
  expect_error(check_progress_args("x.log", NA, "/bin"),
               "must be TRUE or FALSE")
  expect_error(check_progress_args("x.log", c(TRUE, TRUE), "/bin"),
               "must be TRUE or FALSE")
  expect_error(check_progress_args("", FALSE, "/bin"),
               "single file name")
  expect_error(check_progress_args("x.log", FALSE, "remote"),
               "local fit")
  expect_error(check_progress_args("x.log", FALSE, "submit"),
               "local fit")
  expect_error(check_progress_args("no_such_dir/x.log", FALSE, "/bin"),
               "directory of 'progress_file' does not exist")
  expect_error(check_progress_args("definitely_absent.log", TRUE, "/bin"),
               "does not exist")

  ## NULL passes through untouched when not resuming
  expect_null(check_progress_args(NULL, FALSE, "/bin"))
})


test_that("a relative progress_file is resolved against the caller's wd", {

  ## remlf90() setwd()s into tempdir() before running the backend, so a
  ## relative path resolved late would land there and vanish with the session.
  expect_error(pf <- check_progress_args("rel.log", FALSE, "/bin"), NA)
  expect_true(startsWith(pf, normalizePath(getwd(), winslash = "/")))
  expect_identical(basename(pf), "rel.log")
  expect_false(startsWith(pf, normalizePath(tempdir(), winslash = "/")))
})


test_that("resuming is refused under debug", {

  ## The streaming branch is selected on progress_file, not on debug, so a
  ## debug resume would rename the previous log aside, overwrite it, then skip
  ## all parsing and return NULL -- the log consumed for no result.
  expect_error(check_progress_args("x.log", TRUE, "/bin", debug = TRUE),
               "cannot be combined with 'debug = TRUE'")

  ## the combination is only refused together
  expect_error(check_progress_args(NULL, FALSE, "/bin", debug = TRUE), NA)
})


## -- the line index -------------------------------------------------------

test_that("a supplied line index gives identical results", {

  ## extract_block() classifies every line of the log to find the end of a
  ## block. Callers that extract many blocks from one vector can compute that
  ## once; the answer must not change.
  idx <- which(!is_numericlog(x3))
  mk <- grep("^[[:space:]]*new G[[:space:]]*$", x3)[1]

  expect_identical(extract_block(mk + 1L, x3),
                   extract_block(mk + 1L, x3, idx))
  expect_identical(parse_block_strict(mk + 1L, x3),
                   parse_block_strict(mk + 1L, x3, idx))

  a <- last_reml_checkpoint(x3, n_groups = 3L, spd = FALSE)
  expect_equal(attr(a, "round"), 18L)
  for (nm in seq_along(a))
    expect_equal(unname(a[[nm]]),
                 unname(parse.txtmat(extract_block(
                   (grep("Genetic variance|Residual variance", x3) + 1L)[nm], x3))))
})


test_that("the whole log is classified once per walk, not once per block", {

  ## The regression guard for the O(rounds x lines) rescan. Counting sweeps is
  ## deterministic where a wall-clock threshold would flake on a loaded runner.
  ns <- asNamespace("breedR")
  full <- length(y1)
  count <- new.env(parent = emptyenv())

  watch <- function(expr) {
    count$n <- 0L
    ## Name given as a string: trace()/untrace() otherwise resolve it by
    ## non-standard evaluation, which finds an unexported function under
    ## load_all() but not under R CMD check's test_check().
    invisible(capture.output(
      trace("is_numericlog", where = ns, print = FALSE,
            tracer = function() {
              if (length(get("x", envir = parent.frame())) == full)
                count$n <- count$n + 1L
            })))
    on.exit(invisible(capture.output(untrace("is_numericlog", where = ns))),
            add = TRUE)
    force(expr)
    count$n
  }

  ## n_groups it can never satisfy, so all 542 rounds are visited
  swept <- watch(last_reml_checkpoint(y1, n_groups = 9L, spd = FALSE))
  expect_equal(swept, 1L)

  ## for contrast: without a supplied index, five rounds alone cost five sweeps
  anchors <- reml_round_index(y1)
  naive <- watch({
    for (i in utils::tail(seq_along(anchors), 5)) {
      to <- if (i < length(anchors)) anchors[[i + 1L]] - 1L else length(y1)
      reml_round_covariances(y1, anchors[[i]], to)
    }
  })
  expect_gte(naive, 5L)
})



## -- the streamed execution path ------------------------------------------

test_that("remlf90() guards the backend it starts", {

  ## Source-level assertions, deliberately. These behaviours cannot be reached
  ## from a test: remlf90() exposes no handle on the process, there is no
  ## portable way to deliver an interrupt to our own session, and the backend
  ## always ends its output with a newline so the drain never fires in
  ## practice. Without these, deleting any of the three lines leaves the whole
  ## suite green -- which was true of the commit that introduced them.
  src <- paste(deparse(body(remlf90)), collapse = "\n")

  ## the child must be killed on exit, prepended so it runs before the log is
  ## closed, and guarded so a failure cannot skip the handlers after it
  expect_true(grepl("try(px$kill", src, fixed = TRUE))
  expect_true(grepl("after = FALSE", src, fixed = TRUE))

  ## and the buffered tail must still be drained
  expect_true(grepl("tail_out", src, fixed = TRUE))
})


test_that("the tail drain recovers a newline-less line", {

  ## read_output_lines() returns a final line with no newline only once it has
  ## seen EOF. On Windows that usually comes after the loop has already left on
  ## !is_alive(); on Unix it usually comes in time. Whichever path picks the
  ## line up, it must end up in the output exactly once, in order.
  skip_if_not_installed("processx")

  rs <- file.path(R.home("bin"), "Rscript")
  d <- file.path(tempdir(), "drain"); dir.create(d, showWarnings = FALSE)
  f <- file.path(d, "tail.R")
  writeLines("cat(\"alpha\\nbeta\\nno-trailing-newline\")", f)

  px <- processx::process$new(rs, f, stdout = "|")
  out <- character(0)
  repeat {
    px$poll_io(1000)
    l <- px$read_output_lines()
    if (length(l)) out <- c(out, l) else if (!px$is_alive()) break
  }

  tail_out <- px$read_output()
  if (nzchar(tail_out)) out <- c(out, strsplit(tail_out, "\r?\n")[[1]])
  expect_identical(out, c("alpha", "beta", "no-trailing-newline"))
})


test_that("the tail drain splits CRLF the way the line reader does", {

  ## The other thing the drain can pick up is whole lines, when the backend
  ## writes its last block and exits between the line reader coming up empty
  ## and is_alive() being checked. Those arrive with their terminators, and
  ## read_output() does not translate CRLF the way read_output_lines() does --
  ## so splitting on \\n alone would leave a stray CR on every one of them and
  ## make the drained path disagree with the main one.
  crlf <- "alpha\r\nbeta\r\n"

  naive <- strsplit(crlf, "\n", fixed = TRUE)[[1]]
  expect_true(any(grepl("\r", naive)))

  drained <- strsplit(crlf, "\r?\n")[[1]]
  expect_identical(drained, c("alpha", "beta"))
  expect_false(any(grepl("\r", drained)))
})
