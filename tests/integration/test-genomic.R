### Integration test: single-step GBLUP end to end (requires binaries) ###

context("Single-step GBLUP (ssGBLUP)")

## Build a small synthetic genotype set for a subset of globulus animals.
## The SNP file is keyed by the original animal ids; breedR must build the
## PREGSF90 XrefID from its pedigree renumbering for the run to succeed.
set.seed(7)
gen_ids <- globulus$self[1:60]
nsnp <- 200
Gmat <- matrix(sample(0:2, length(gen_ids) * nsnp, replace = TRUE,
                      prob = c(0.25, 0.5, 0.25)),
               nrow = length(gen_ids))
snp_file <- file.path(tempdir(), "test_ssgblup_snp.txt")
write_snp_file(Gmat, ids = gen_ids, file = snp_file)
file.remove(list.files(dirname(snp_file), pattern = "_XrefID$",
                       full.names = TRUE))

res.gen <- suppressMessages(
  remlf90(fixed = phe_X ~ gg,
          genetic = list(model = 'add_animal',
                         pedigree = globulus[, 1:3], id = 'self'),
          genomic = list(snp_file = snp_file, verify_parentage = 0L),
          data = globulus))

test_that("remlf90(genomic=) fits a single-step GBLUP model", {
  expect_s3_class(res.gen, "remlf90")
  expect_true(is.finite(as.numeric(logLik(res.gen))))
  expect_true("gg" %in% names(fixef(res.gen)))
})

test_that("genomic QC results are attached and consistent", {
  expect_false(is.null(res.gen$genomic))
  qc <- res.gen$genomic
  expect_equal(qc$n_snp_total, nsnp)
  expect_true(qc$n_snp_passed >= 0 && qc$n_snp_passed <= nsnp)
  expect_equal(qc$n_snp_excluded, qc$n_snp_total - qc$n_snp_passed)
  expect_true(all(c("snp", "frequency") %in% names(qc$freq)))
})

test_that("a second fit reuses the staged genotype file", {

  ## Each fit works in its own directory, so the genotype file the backend
  ## needs to find there was copied once per fit and never reclaimed -- a
  ## cross-validation loop over a large SNP file would fill the temp volume.
  ## It is staged once per session now, and hard linked into each fit.
  res.gen2 <- suppressMessages(
    remlf90(fixed = phe_X ~ 1,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = snp_file, verify_parentage = 0L),
            data = globulus))

  expect_false(identical(res.gen$reml$dir, res.gen2$reml$dir))

  ## both fits found the genotypes ...
  for (r in list(res.gen, res.gen2))
    expect_true(file.exists(file.path(r$reml$dir, basename(snp_file))))

  ## ... and there is still exactly one staged copy behind them
  cache <- breedR_input_cache(snp_file)
  expect_length(list.files(cache), 1L)

  ## on a filesystem with hard links the two fits share one inode, so the
  ## second fit costs no extra bytes; where they are unavailable stage_input()
  ## falls back to copying and this is merely equality of content
  expect_identical(
    readLines(file.path(res.gen2$reml$dir, basename(snp_file))),
    readLines(file.path(cache, basename(snp_file))))

  clean_workdir(res.gen2)
  expect_false(dir.exists(res.gen2$reml$dir))
})

test_that("a genotyped animal absent from the pedigree is an error", {
  bad_file <- file.path(tempdir(), "test_ssgblup_bad.txt")
  write_snp_file(Gmat[1:3, , drop = FALSE],
                 ids = c(9999999, gen_ids[2:3]), file = bad_file)
  file.remove(list.files(dirname(bad_file), pattern = "_XrefID$",
                         full.names = TRUE))
  expect_error(
    suppressMessages(
      remlf90(fixed = phe_X ~ gg,
              genetic = list(model = 'add_animal',
                             pedigree = globulus[, 1:3], id = 'self'),
              genomic = list(snp_file = bad_file, verify_parentage = 0L),
              data = globulus)),
    "not found in the pedigree")
})


context("Standalone genotype QC (qcf90)")

test_that("qcf90 runs and returns clean-marker results", {
  qc_snp <- file.path(tempdir(), "test_qcf90_snp.txt")
  write_snp_file(Gmat, ids = gen_ids, file = qc_snp)
  qc <- qcf90(snp_file = qc_snp)
  expect_type(qc, "list")
  expect_false(is.null(qc$clean_snp))    # a cleaned SNP file was produced
  expect_false(is.null(qc$log))
})


context("Genome-wide association (ssGWAS via postgsf90)")

test_that("postgsf90 backsolves SNP effects on a save_ginverse fit", {
  gw_snp <- file.path(tempdir(), "test_gwas_snp.txt")
  write_snp_file(Gmat, ids = gen_ids, file = gw_snp)
  file.remove(list.files(dirname(gw_snp), pattern = "_XrefID$",
                         full.names = TRUE))
  res.gwas <- suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = gw_snp, verify_parentage = 0L,
                           save_ginverse = TRUE),
            data = globulus))
  gwas <- postgsf90(res.gwas, manhattan_plot = TRUE)
  expect_equal(nrow(gwas$snp_sol), nsnp)     # one SNP solution per marker
  expect_false(is.null(gwas$manhattan))
  expect_equal(nrow(gwas$manhattan), nsnp)
})

test_that("postgsf90 refuses a fit made without save_ginverse", {
  ng_snp <- file.path(tempdir(), "test_nogwas_snp.txt")
  write_snp_file(Gmat, ids = gen_ids, file = ng_snp)
  file.remove(list.files(dirname(ng_snp), pattern = "_XrefID$",
                         full.names = TRUE))
  res.ng <- suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = ng_snp, verify_parentage = 0L),
            data = globulus))
  expect_error(postgsf90(res.ng), "save_ginverse")
})


context("Genomic prediction of DGV (predf90)")

test_that("predf90 predicts DGV from postgsf90 SNP effects", {
  # In practice predf90 predicts for NEW animals; here we reuse the genotyped
  # set as an end-to-end smoke test (predf90 needs postgsf90's snp_pred).
  gp_snp <- file.path(tempdir(), "test_predf90_snp.txt")
  write_snp_file(Gmat, ids = gen_ids, file = gp_snp)
  file.remove(list.files(dirname(gp_snp), pattern = "_XrefID$",
                         full.names = TRUE))
  res.gp <- suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = gp_snp, verify_parentage = 0L,
                           save_ginverse = TRUE),
            data = globulus))
  gwas <- postgsf90(res.gp, manhattan_plot = TRUE)

  ## postgsf90() reports where it wrote the SNP effects; predf90() has to be
  ## pointed at that, since each fit works in its own directory.
  expect_identical(gwas$dir, res.gp$reml$dir)

  pred <- predf90(snp_file = gp_snp, dir = gwas$dir)
  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(Gmat))          # one DGV per genotyped animal
  expect_true(all(c("id", "dgv") %in% names(pred)))
  expect_true(all(is.finite(pred$dgv)))
})

test_that("predf90 on a validation set leaves the training genotypes alone", {

  ## The realistic call: predict for animals the model was not fitted on, from
  ## a genotype file the user happens to have named the same thing. The staged
  ## file in a fit's directory is a hard link into the session input cache, so
  ## copying the validation file over it wrote through the link and replaced
  ## the training genotypes in the cache -- and so in every other fit directory
  ## of the session -- with no error anywhere.
  train_dir <- breedR_workdir('train_')
  valid_dir <- breedR_workdir('valid_')
  sibling   <- breedR_workdir('sibling_')
  on.exit(unlink(c(train_dir, valid_dir, sibling), recursive = TRUE),
          add = TRUE)

  ## Same name, same dimensions, hence the same byte count: the two files are
  ## distinguishable only by their contents.
  train_geno <- file.path(train_dir, "genotypes.txt")
  valid_geno <- file.path(valid_dir, "genotypes.txt")
  write_snp_file(Gmat, ids = gen_ids, file = train_geno)
  Vmat <- matrix(sample(0:2, length(gen_ids) * nsnp, replace = TRUE,
                        prob = c(0.25, 0.5, 0.25)),
                 nrow = length(gen_ids))
  write_snp_file(Vmat, ids = gen_ids, file = valid_geno)
  file.remove(list.files(c(train_dir, valid_dir), pattern = "_XrefID$",
                         full.names = TRUE))
  expect_identical(file.info(train_geno)$size, file.info(valid_geno)$size)

  train_lines <- readLines(train_geno)

  res.vp <- suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = train_geno, verify_parentage = 0L,
                           save_ginverse = TRUE),
            data = globulus))
  gwas <- postgsf90(res.vp)

  ## A second directory staged from the same source, standing in for another
  ## fit of the session. Resolve the cache directory now: asking for it after
  ## the fact would re-stage a damaged copy and hide the very thing under test.
  stage_input(train_geno, sibling)
  cached <- file.path(breedR_input_cache(train_geno), basename(train_geno))

  pred <- predf90(snp_file = valid_geno, dir = gwas$dir)
  expect_s3_class(pred, "data.frame")

  ## The fit's own directory now holds the validation genotypes, as asked ...
  expect_identical(readLines(file.path(gwas$dir, "genotypes.txt")),
                   readLines(valid_geno))

  ## ... and nothing reachable through the link went with it.
  expect_identical(readLines(cached), train_lines)
  expect_identical(readLines(file.path(sibling, "genotypes.txt")), train_lines)
})

test_that("predf90 errors without a prior postgsf90 run", {
  lone_snp <- file.path(tempdir(), "test_predf90_lone.txt")
  write_snp_file(Gmat, ids = gen_ids, file = lone_snp)
  expect_error(
    predf90(snp_file = lone_snp, dir = file.path(tempdir(), "no_postgs")),
    "snp_pred")
})
