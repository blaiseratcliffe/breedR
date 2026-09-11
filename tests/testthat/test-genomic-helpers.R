## Tests for genomic helper functions (R/genomic.R)

context("Genomic helper functions")

## -- check_genomic() --

test_that("check_genomic requires snp_file", {
  expect_error(check_genomic(list()), "snp_file.*required")
})

test_that("check_genomic rejects non-existent snp_file", {
  expect_error(check_genomic(list(snp_file = "nonexistent.txt")),
               "not found")
})

test_that("check_genomic rejects non-existent map_file", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), map_file = "no.txt")),
    "not found")
})

test_that("check_genomic sets correct defaults", {
  result <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR")))
  expect_equal(result$whichG, 1L)
  expect_equal(result$tunedG, 2L)
  expect_equal(result$AlphaBeta, c(0.95, 0.05))
  expect_equal(result$minfreq, 0.05)
  expect_equal(result$callrate, 0.90)
  expect_equal(result$callrateAnim, 0.90)
  expect_equal(result$verify_parentage, 3L)
  expect_false(result$saveG)
  expect_false(result$saveA22)
})

test_that("check_genomic rejects invalid whichG", {
  expect_error(check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), whichG = 5)),
               "whichG.*must be 1, 2, or 3")
})

test_that("check_genomic rejects invalid tunedG", {
  expect_error(check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), tunedG = 7)),
               "tunedG")
})

test_that("check_genomic rejects invalid AlphaBeta", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), AlphaBeta = 0.95)),
    "length 2")
})

test_that("check_genomic rejects invalid callrate", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), callrate = 1.5)),
    "between 0 and 1")
})

test_that("check_genomic rejects invalid verify_parentage", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), verify_parentage = 5)),
    "verify_parentage")
})

test_that("check_genomic accepts valid extra_options", {
  result <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                extra_options = c("no_quality_control")))
  expect_equal(result$extra_options, "no_quality_control")
})

test_that("check_genomic rejects non-character extra_options", {
  expect_error(
    check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"), extra_options = 42)),
    "character vector")
})

## -- build_genomic_options() --

test_that("build_genomic_options generates correct options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR")))
  opts <- build_genomic_options(genomic)

  expect_true(any(grepl("^SNP_file DESCRIPTION$", opts)))
  expect_true(any(grepl("^whichG 1$", opts)))
  expect_true(any(grepl("^tunedG 2$", opts)))
  expect_true(any(grepl("^AlphaBeta 0\\.95 0\\.05$", opts)))
  expect_true(any(grepl("^minfreq 0\\.05$", opts)))
  expect_true(any(grepl("^callrate 0\\.9$", opts)))
  expect_true(any(grepl("^verify_parentage 3$", opts)))
})

test_that("build_genomic_options includes map_file when provided", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 map_file = system.file("NAMESPACE", package = "breedR")))
  opts <- build_genomic_options(genomic)
  expect_true(any(grepl("^map_file", opts)))
})

test_that("build_genomic_options includes save options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 saveG = TRUE, saveA22 = TRUE))
  opts <- build_genomic_options(genomic)
  expect_true("saveG" %in% opts)
  expect_true("saveA22" %in% opts)
})

test_that("build_genomic_options includes extra_options", {
  genomic <- check_genomic(list(snp_file = system.file("DESCRIPTION", package = "breedR"),
                                 extra_options = c("no_quality_control",
                                                   "thrStopCorAG 0.0")))
  opts <- build_genomic_options(genomic)
  expect_true("no_quality_control" %in% opts)
  expect_true("thrStopCorAG 0.0" %in% opts)
})


## -- stage_input() --

## Each fit works in its own directory, so the genotype file the backend needs
## to find there used to be copied once per fit and never reclaimed. It is now
## copied into a session cache once and linked from there.

make_src <- function(dir, name = 'geno.txt', text = 'a b c') {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, name)
  writeLines(text, f)
  f
}

test_that("stage_input() puts the file in the working directory", {

  src <- make_src(breedR_workdir('src_'))
  wd  <- breedR_workdir()
  on.exit(unlink(c(dirname(src), wd), recursive = TRUE), add = TRUE)

  stage_input(src, wd)

  dest <- file.path(wd, basename(src))
  expect_true(file.exists(dest))
  expect_identical(readLines(dest), readLines(src))
})


test_that("stage_input() keeps one cached copy across working directories", {

  src <- make_src(breedR_workdir('src_'))
  w1  <- breedR_workdir()
  w2  <- breedR_workdir()
  on.exit(unlink(c(dirname(src), w1, w2), recursive = TRUE), add = TRUE)

  stage_input(src, w1)
  c1 <- breedR_input_cache(src)
  stage_input(src, w2)
  c2 <- breedR_input_cache(src)

  ## both fits can read it ...
  expect_true(file.exists(file.path(w1, basename(src))))
  expect_true(file.exists(file.path(w2, basename(src))))

  ## ... off a single cached copy. This is the whole fix: without it a session
  ## of N genomic fits held N copies of the genotype file.
  expect_identical(c1, c2)
  expect_length(list.files(c1), 1L)
})


test_that("stage_input() re-stages a source that changed on disk", {

  ## A cached copy is reused while the source's size and mtime are unchanged,
  ## so this is the case that would serve stale genotypes if the entry were
  ## keyed on the path alone.
  src <- make_src(breedR_workdir('src_'), text = 'first')
  wd  <- breedR_workdir()
  on.exit(unlink(c(dirname(src), wd), recursive = TRUE), add = TRUE)

  stage_input(src, wd)
  expect_identical(readLines(file.path(wd, basename(src))), 'first')

  writeLines(c('second', 'and longer'), src)
  Sys.setFileTime(src, Sys.time() + 5)
  stage_input(src, wd)

  expect_identical(readLines(file.path(wd, basename(src))),
                   c('second', 'and longer'))
})


test_that("stage_input() leaves a file that is already in place alone", {

  wd  <- breedR_workdir()
  src <- make_src(wd)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)

  expect_silent(stage_input(src, wd))
  expect_true(file.exists(src))
  expect_identical(readLines(src), 'a b c')
})


test_that("stage_input() re-stages when the cached copy itself changed", {

  ## The staged file is hard linked into every working directory, so a backend
  ## that rewrote its input in place would reach back through the link and
  ## alter the cache. Comparing only against the source would never see it.
  src <- make_src(breedR_workdir('src_'), text = 'genotypes')
  wd  <- breedR_workdir()
  on.exit(unlink(c(dirname(src), wd), recursive = TRUE), add = TRUE)

  stage_input(src, wd)
  cached <- file.path(breedR_input_cache(src), basename(src))

  writeLines('clobbered by the backend', cached)
  stage_input(src, wd)

  expect_identical(readLines(file.path(wd, basename(src))), 'genotypes')
})


test_that("stage_input() re-stages when the cached copy changed but kept its size", {

  ## The case above changes the length too, so a size comparison alone catches
  ## it. The rewrite that actually happens does not: two genotype files for the
  ## same markers are the same length, and a backend rewriting fixed width
  ## records in place keeps the length by construction. Only the copy's own
  ## mtime separates this from an untouched cache.
  src <- make_src(breedR_workdir('src_'), text = 'genotypes')
  wd  <- breedR_workdir()
  on.exit(unlink(c(dirname(src), wd), recursive = TRUE), add = TRUE)

  stage_input(src, wd)
  cached <- file.path(breedR_input_cache(src), basename(src))
  before <- file.info(cached)$size

  writeLines('clobbered', cached)          # same nine characters
  expect_identical(file.info(cached)$size, before)

  ## Force the mtime apart rather than trusting the filesystem to resolve two
  ## writes a moment apart, as the source-side test above does.
  Sys.setFileTime(cached, Sys.time() + 5)

  stage_input(src, wd)

  expect_identical(readLines(file.path(wd, basename(src))), 'genotypes')
})


test_that("stage_input() does not write through the link into sibling directories", {

  ## Staging a *different* source that happens to share a basename -- a
  ## validation genotype file against a directory holding the training one,
  ## which is what predf90() does -- must replace the link, not write through
  ## it. Writing through reaches the cache entry and so every other working
  ## directory staged from the same source.
  train <- make_src(breedR_workdir('train_'), text = 'training genotypes')
  valid <- make_src(breedR_workdir('valid_'), text = 'validate genotypes')
  w1    <- breedR_workdir()
  w2    <- breedR_workdir()
  on.exit(unlink(c(dirname(train), dirname(valid), w1, w2), recursive = TRUE),
          add = TRUE)

  expect_identical(basename(train), basename(valid))

  stage_input(train, w1)
  stage_input(train, w2)

  stage_input(valid, w1)

  ## The directory asked for the swap sees it ...
  expect_identical(readLines(file.path(w1, basename(valid))),
                   'validate genotypes')

  ## ... and nothing else does.
  expect_identical(readLines(file.path(w2, basename(train))),
                   'training genotypes')
  expect_identical(readLines(file.path(breedR_input_cache(train),
                                       basename(train))),
                   'training genotypes')
})


## -- write_xref_from_pedigree() --

test_that("write_xref_from_pedigree() finds genotyped ids of 1e5 and above (#43)", {

  ## The SNP file holds ids as text, so they are matched against the pedigree
  ## labels as text. A double code 100000 used to be labelled "1e+05".
  ped <- build_pedigree(1:3, data = data.frame(self = as.numeric(seq_len(1e5)),
                                               dad = 0, mum = 0))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')
  writeLines(c('100000 0120', '5 1111'), snp)

  write_xref_from_pedigree(snp, ped, wd)

  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   c('100000 100000', '5 5'))
})


## -- postgsf90() input checks --

## These run without binaries: each stops before anything is executed.

test_that("postgsf90() rejects a model it cannot use", {

  expect_error(postgsf90(list()), "must be a fitted remlf90 object")

  no_genomic <- structure(list(), class = c('breedR', 'remlf90'))
  expect_error(postgsf90(no_genomic), "fitted with genomic")

  no_ginv <- structure(list(genomic = list(save_ginverse = FALSE)),
                       class = c('breedR', 'remlf90'))
  expect_error(postgsf90(no_ginv), "not fitted for GWAS")
})


test_that("postgsf90() says so when the model records no working directory", {

  ## A local fit always records one, so an object without the field never had
  ## it: an older version of the package, or a job recovered through
  ## breedR.qget(), which rebuilds the object through parse_results() and skips
  ## the assignment in remlf90(). It used to fall back to bare tempdir(), which
  ## holds nobody's solutions, so the report was that the fit was missing from
  ## a directory it had never been in.
  stale <- structure(list(genomic = list(save_ginverse = TRUE)),
                     class = c('breedR', 'remlf90'))

  msg <- tryCatch(postgsf90(stale), error = conditionMessage)

  expect_match(msg, "no recorded working directory")
  expect_match(msg, "Refit it in this session")
  expect_match(msg, "recovered from a submitted job")

  ## and it does not go looking in the session's temporary files first
  expect_false(grepl("Solutions file not found", msg, fixed = TRUE))

  ## It used to blame a saved file, which is the one case that cannot get
  ## here: saving keeps the field. See the test below.
  expect_false(grepl("saved file", msg, fixed = TRUE))
})


test_that("postgsf90() tells a reloaded model apart from one that never had a directory", {

  ## saveRDS() keeps reml$dir -- what a reload loses is the directory it names,
  ## not the field. So this is the case that reaches the *second* check, and
  ## the message above, which blames an old version of the package, would be
  ## the wrong thing to say about it.
  gone <- file.path(tempdir(), "a_session_that_ended")
  reloaded <- structure(list(genomic = list(save_ginverse = TRUE),
                             reml = list(dir = gone)),
                        class = c('breedR', 'remlf90'))

  msg <- tryCatch(postgsf90(reloaded), error = conditionMessage)

  expect_match(msg, "Solutions file not found")
  expect_match(msg, gone, fixed = TRUE)
  expect_match(msg, "same R session")

  ## The path named above is from a session that has ended, so say so rather
  ## than leaving a stale absolute path to read as a bug in the lookup.
  expect_match(msg, "went away with it")

  ## The directory in the message is named, so it should not also be described
  ## as never having been recorded.
  expect_false(grepl("no recorded working directory", msg, fixed = TRUE))
})


test_that("predf90() requires the directory postgsf90() wrote to", {

  ## Nothing writes snp_pred into bare tempdir() now that each fit works in its
  ## own subdirectory, so the old default could only ever end at "snp_pred file
  ## not found". Better to ask for the one value that can be right.
  snp <- make_src(breedR_workdir('predsrc_'), name = 'geno.txt')
  on.exit(unlink(dirname(snp), recursive = TRUE), add = TRUE)

  expect_error(predf90(snp_file = snp), 'argument "dir" is missing')
})
