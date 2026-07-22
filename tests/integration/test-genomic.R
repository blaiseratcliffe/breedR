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
