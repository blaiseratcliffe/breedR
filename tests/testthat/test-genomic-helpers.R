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

test_that("write_xref_from_pedigree() finds ids written by write_snp_file() from double ids of 1e5 and above (#54)", {

  ## write_snp_file() used to write ids = 1e5 (a double) as "1e+05"
  ## (as.character() on a double), and no genotyped animal matched the
  ## pedigree integer-formatted labels. Same quirk as #43, in the SNP
  ## writer rather than build_pedigree().
  ped <- build_pedigree(1:3, data = data.frame(self = as.numeric(seq_len(1e5)),
                                               dad = 0, mum = 0))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')
  write_snp_file(matrix(c(0L, 1L, 2L), nrow = 1), ids = 1e5, file = snp)

  write_xref_from_pedigree(snp, ped, wd)

  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   '100000 100000')
})

test_that("write_xref_from_pedigree() translates ids of a recoded pedigree (#49)", {

  ## Animal 1 precedes its parents 2 and 3, so the pedigree is recoded with
  ## map 3 1 2. The SNP file holds the original ids; the XrefID must give the
  ## recoded code that breedR writes in the pedigree and data files.
  ped <- suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = 1:3, dad = c(2L, 0L, 0L),
                                          mum = c(3L, 0L, 0L))))
  expect_identical(attr(ped, 'map'), c(3L, 1L, 2L))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')
  writeLines(c('1 0120', '2 1111', '3 2222'), snp)

  write_xref_from_pedigree(snp, ped, wd)

  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   c('3 1', '1 2', '2 3'))
})

test_that("write_xref_from_pedigree() rejects a recoded code given as an id (#49)", {

  ## Codes with gaps are recoded 10, 20, 30 -> 1, 2, 3. An id of 2 is not an
  ## animal of this pedigree, even though 2 is one of its recoded codes.
  ped <- suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = c(10L, 20L, 30L),
                                          dad = 0L, mum = 0L)))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')

  writeLines(c('30 0120', '10 1111'), snp)
  write_xref_from_pedigree(snp, ped, wd)
  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   c('3 30', '1 10'))

  writeLines(c('30 0120', '2 1111'), snp)
  expect_error(write_xref_from_pedigree(snp, ped, wd),
               "not found in the pedigree: 2$")
})

test_that("write_xref_from_pedigree() finds ids of 1e5 and above in a recoded pedigree (#43, #49)", {

  ## Double codes with gaps are recoded. The ids translated back through the
  ## map must print as "100000", not "1e+05", or no genotyped animal matches.
  ped <- suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = c(1e5, 2.5e5, 300001),
                                          dad = 0, mum = 0)))
  expect_false(is.null(attr(ped, 'map')))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')
  writeLines(c('300001 0120', '100000 1111'), snp)

  write_xref_from_pedigree(snp, ped, wd)

  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   c('3 300001', '1 100000'))
})

test_that("gibbsf90(genomic=) derives the XrefID from its pedigree (#55)", {

  ## 1..3 come before their parents 4..6, so the pedigree is recoded with
  ## map 3 4 6 5 1 2; each code is also another animal's id
  ped <- data.frame(self = 1:6, dad = c(5L, 5L, 5L, 0L, 0L, 0L),
                    mum = c(6L, 6L, 4L, 0L, 0L, 0L))
  dat <- data.frame(ped[1:4, ], y = c(1.2, 0.3, 2.1, 1.7))
  expect_identical(
    attr(suppressWarnings(build_pedigree(1:3, data = ped)), 'map'),
    c(3L, 4L, 6L, 5L, 1L, 2L))
  wd <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  bin <- file.path(wd, 'bin')
  dir.create(bin)
  file.create(file.path(bin, grep('gibbsf90+', genomic_program_files(),
                                  value = TRUE, fixed = TRUE)))
  snp <- file.path(wd, 'geno.txt')
  write_snp_file(matrix(c(0L,1L,2L,1L, 2L,1L,0L,1L, 1L,1L,1L,0L), 3),
                 ids = c(1L, 3L, 5L), file = snp)
  captured <- new.env()
  local_mocked_bindings(
    check_genomic_programs = function(...) TRUE,
    run_pregsf90 = function(dir, bin_path) {
      captured$dir  <- dir
      captured$xref <- readLines(file.path(dir, 'geno.txt_XrefID'))
      stop('preGSf90 reached')
    }, .package = 'breedR')
  expect_error(suppressWarnings(
    gibbsf90(y ~ 1,
             genetic = list(model = 'add_animal', pedigree = ped, id = 'self'),
             genomic = list(snp_file = snp), data = dat, breedR.bin = bin,
             n_samples = 10L)),
    'preGSf90 reached')
  if (!is.null(captured$dir)) unlink(captured$dir, recursive = TRUE)
  expect_identical(captured$xref, c('3 1', '6 3', '1 5'))
})

test_that("write_xref_from_pedigree() translates ids of a pedigree coded from 2 (#50)", {

  ## Codes 2..4 are recoded to 1..3, so the map starts with NA. The genotyped
  ## animals must still get the recoded code of their own original id.
  ped <- suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = 2:4, dad = c(0L, 0L, 2L),
                                          mum = c(0L, 0L, 3L))))
  expect_identical(attr(ped, 'map'), c(NA, 1L, 2L, 3L))
  wd  <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  snp <- file.path(wd, 'geno.txt')
  writeLines(c('4 0120', '2 1111'), snp)

  write_xref_from_pedigree(snp, ped, wd)

  expect_identical(readLines(file.path(wd, 'geno.txt_XrefID')),
                   c('3 4', '1 2'))
})


## -- a supplied <snp_file>_XrefID --

## preGSf90 pairs XrefID row i with SNP-file row i and takes the animal's code
## from column 1; it never reads column 2. Runs the genomic pipeline on the
## SNP file '1 / 2 / 3' with the given XrefID lines next to it, preGSf90
## mocked. Returns the pipeline's error, the bytes of the XrefID staged for
## preGSf90 (NULL if it was never reached) and the bytes supplied.
run_with_supplied_xref <- function(xref_lines, pedigree) {
  src <- breedR_workdir()
  wd  <- breedR_workdir()
  on.exit(unlink(c(src, wd), recursive = TRUE), add = TRUE)
  snp <- file.path(src, 'geno.txt')
  writeLines(c('1 0120', '2 1111', '3 2222'), snp)
  xref <- paste0(snp, '_XrefID')
  writeLines(xref_lines, xref)
  staged <- NULL
  local_mocked_bindings(
    write.progsf90 = function(...) invisible(NULL),
    check_genomic_programs = function(...) TRUE,
    run_pregsf90 = function(dir, bin_path) {
      f <- file.path(dir, 'geno.txt_XrefID')
      staged <<- readBin(f, 'raw', file.size(f))
      stop('preGSf90 reached')
    }, .package = 'breedR')
  err <- tryCatch(
    run_pregsf90_pipeline(list(snp_file = snp), character(0),
                          list(parameter = list(options = NULL)), wd, wd,
                          pedigree = pedigree),
    error = conditionMessage)
  list(error = err, staged = staged,
       supplied = readBin(xref, 'raw', file.size(xref)))
}

## Animal 1 precedes its parents 2 and 3, so the pedigree is recoded with map
## 3 1 2: the XrefID breedR derives for the SNP file above is '3 1 / 1 2 / 2 3'.
recoded_ped_1to3 <- function()
  suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = 1:3, dad = c(2L, 0L, 0L),
                                          mum = c(3L, 0L, 0L))))

test_that("a supplied XrefID with its rows out of order is refused (#53)", {
  ped <- recoded_ped_1to3()
  expect_identical(attr(ped, 'map'), c(3L, 1L, 2L))

  res <- run_with_supplied_xref(c('1 2', '3 1', '2 3'), ped)
  expect_match(res$error, "Row 1 of 'geno.txt_XrefID'", fixed = TRUE)
  expect_null(res$staged)
})

test_that("a supplied XrefID with codes other than breedR's is refused (#53)", {
  ped <- recoded_ped_1to3()

  ## original ids in column 1, the #49 mistake on a recoded pedigree
  res <- run_with_supplied_xref(c('1 1', '2 2', '3 3'), ped)
  expect_match(res$error,
               "Row 1 of 'geno.txt_XrefID' codes animal '1' as 1, but breedR codes it 3",
               fixed = TRUE)
  expect_null(res$staged)

  ## the right ids in column 2, the codes in reverse order
  res <- run_with_supplied_xref(c('2 1', '1 2', '3 3'), ped)
  expect_match(res$error,
               "Row 1 of 'geno.txt_XrefID' codes animal '1' as 2, but breedR codes it 3",
               fixed = TRUE)
  expect_null(res$staged)
})

test_that("a supplied XrefID with fewer rows than the SNP file is refused (#53)", {
  res <- run_with_supplied_xref(c('3 1', '1 2'), recoded_ped_1to3())
  expect_match(res$error, "'geno.txt_XrefID' has 2 rows", fixed = TRUE)
  expect_match(res$error, "genotype file 3", fixed = TRUE)
  expect_null(res$staged)
})

test_that("a supplied XrefID out of order is refused without a pedigree (#53)", {
  res <- run_with_supplied_xref(c('1 2', '3 1', '2 3'), NULL)
  expect_match(res$error, "Row 1 of 'geno.txt_XrefID'", fixed = TRUE)
  expect_null(res$staged)
})

test_that("a correct supplied XrefID reaches preGSf90 byte for byte (#53)", {
  ped <- recoded_ped_1to3()

  res <- run_with_supplied_xref(c('3 1', '1 2', '2 3'), ped)
  expect_identical(res$error, 'preGSf90 reached')
  expect_identical(res$staged, res$supplied)

  res <- run_with_supplied_xref(c('  3   1', '1\t2 ', '2 3'), ped)
  expect_identical(res$error, 'preGSf90 reached')
  expect_identical(res$staged, res$supplied)
})


## -- parse_pregsf90_qc() --

## preGSf90 reports animals by the codes of the pedigree breedR wrote for it,
## which are breedR's internal codes. The reports below are laid out as
## preGSf90 writes them.
write_qc_reports <- function(dir) {
  writeLines(
    " Genotyped Animal with low call rate REMOVED         3    0.2500",
    file.path(dir, 'Gen_call_rate'))
  writeLines(c(
    " Quality Control - Check Parent-Progeny Mendelian  conflicts",
    "",
    "    Total animals: 3 - Genotyped animals: 3 - Effective: 2",
    "",
    "   Animal - Sire Conflict         3         1        81    0.4050",
    "   Animal - Dam Conflict         3         2        24    0.1200",
    "",
    "   Number of Parent-Progeny Mendelian Conflicts: 2",
    " ",
    "       #_Gen   Renf90_Id      #_sire  tot_#_sire       #_dam   tot_#_dam       #_ind   tot_#_ind",
    "           2           3           0           0           0           0           2           2",
    "           1           1           1           1           0           0           0           0"),
    file.path(dir, 'Gen_conflicts'))
}

test_that("parse_pregsf90_qc() names animals of a recoded pedigree by their ids (#58)", {

  ## Animal 1 precedes its parents 2 and 3, so the pedigree is recoded with
  ## map 3 1 2: animal 1 is code 3, animal 2 code 1 and animal 3 code 2. Every
  ## code is also another animal's id, so an untranslated code names the
  ## wrong animal.
  ped <- suppressWarnings(
    build_pedigree(1:3, data = data.frame(self = 1:3, dad = c(2L, 0L, 0L),
                                          mum = c(3L, 0L, 0L))))
  expect_identical(attr(ped, 'map'), c(3L, 1L, 2L))
  wd <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  write_qc_reports(wd)

  qc <- parse_pregsf90_qc(wd, pedigree = ped)

  ## code 3 is animal 1
  expect_identical(qc$excluded_animals[[8]], 1L)
  expect_identical(qc$excluded_animals[[9]], 0.25)
  expect_identical(qc$n_animals_excluded, 1L)

  ## Progeny 1 with its sire 2 and its dam 3. In the table only the second
  ## column is a pedigree code: the first is the row of the SNP file, and it
  ## must stay 2 and 1. Every other line is left as preGSf90 wrote it.
  expect_identical(qc$conflicts, c(
    " Quality Control - Check Parent-Progeny Mendelian  conflicts",
    "",
    "    Total animals: 3 - Genotyped animals: 3 - Effective: 2",
    "",
    "   Animal - Sire Conflict         1         2        81    0.4050",
    "   Animal - Dam Conflict         1         3        24    0.1200",
    "",
    "   Number of Parent-Progeny Mendelian Conflicts: 2",
    " ",
    "       #_Gen   Renf90_Id      #_sire  tot_#_sire       #_dam   tot_#_dam       #_ind   tot_#_ind",
    "           2           1           0           0           0           0           2           2",
    "           1           2           1           1           0           0           0           0"))
})

test_that("parse_pregsf90_qc() leaves the report of a pedigree that was not recoded alone (#58)", {

  p0 <- build_pedigree(1:3, data = data.frame(self = 1:4, dad = c(0, 0, 1, 1),
                                              mum = c(0, 0, 2, 2)))
  expect_null(attr(p0, 'map'))
  wd <- breedR_workdir()
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  write_qc_reports(wd)

  qc <- parse_pregsf90_qc(wd, pedigree = p0)
  expect_identical(qc, parse_pregsf90_qc(wd))
  expect_identical(qc$excluded_animals[[8]], 3L)
  expect_identical(qc$conflicts, readLines(file.path(wd, 'Gen_conflicts')))
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
