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

test_that("relabelling the animals so the pedigree is recoded changes nothing (#49)", {

  ## Reversing the codes puts offspring before their parents, so breedR
  ## recodes the pedigree. The same animals, phenotypes and genotypes under
  ## other names must give the same fit; with the genotypes attached to the
  ## wrong animals they did not.
  N   <- max(globulus[, c('self', 'dad', 'mum')])
  flip <- function(x) ifelse(x == 0, 0L, as.integer(N + 1L - x))
  glob_rev <- globulus
  glob_rev[, c('self', 'dad', 'mum')] <-
    lapply(globulus[, c('self', 'dad', 'mum')], flip)

  rev_dir  <- breedR_workdir('recoded_')
  on.exit(unlink(rev_dir, recursive = TRUE), add = TRUE)
  rev_file <- file.path(rev_dir, "test_ssgblup_recoded.txt")
  write_snp_file(Gmat, ids = flip(gen_ids), file = rev_file)

  res.rev <- suppressWarnings(suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = glob_rev[, 1:3], id = 'self'),
            genomic = list(snp_file = rev_file, verify_parentage = 0L),
            data = glob_rev)))
  map <- attr(get_pedigree(res.rev), 'map')
  expect_false(is.null(map))

  ## Relabelling moves logLik by a relative 1e-14 or so, with the Windows and
  ## the Linux binaries. With the genotypes on the wrong animals it moved by a
  ## relative 1.8e-6.
  expect_equal(as.numeric(logLik(res.rev)), as.numeric(logLik(res.gen)),
               tolerance = 1e-9)
  expect_equal(res.rev$var, res.gen$var, tolerance = 1e-6)
  expect_equal(fitted(res.rev), fitted(res.gen), tolerance = 1e-6)

  ## Each animal's breeding value, matched by its globulus id. The recoded fit
  ## names them by the ids it was given, which flip() takes back.
  bv_rev <- ranef(res.rev)$genetic
  bv_gen <- ranef(res.gen)$genetic
  orig   <- flip(as.integer(names(bv_rev)))
  expect_equal(as.numeric(bv_rev)[match(as.integer(names(bv_gen)), orig)],
               as.numeric(bv_gen), tolerance = 1e-6)
})

test_that("codes that start above 1 are recoded and change nothing (#50)", {

  ## Adding 100 to every code keeps the pedigree sorted and consecutive, but it
  ## no longer starts at 1. It used to go unrecoded, the codes were taken as
  ## positions, and the fit stopped.
  up <- function(x) ifelse(x == 0, 0L, as.integer(x + 100L))
  glob_up <- globulus
  glob_up[, c('self', 'dad', 'mum')] <-
    lapply(globulus[, c('self', 'dad', 'mum')], up)

  up_dir  <- breedR_workdir('shifted_')
  on.exit(unlink(up_dir, recursive = TRUE), add = TRUE)
  up_file <- file.path(up_dir, "test_ssgblup_shifted.txt")
  write_snp_file(Gmat, ids = up(gen_ids), file = up_file)

  res.up <- suppressWarnings(suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = glob_up[, 1:3], id = 'self'),
            genomic = list(snp_file = up_file, verify_parentage = 0L),
            data = glob_up)))
  map <- attr(get_pedigree(res.up), 'map')
  expect_false(is.null(map))

  expect_equal(as.numeric(logLik(res.up)), as.numeric(logLik(res.gen)),
               tolerance = 1e-9)
  expect_equal(res.up$var, res.gen$var, tolerance = 1e-6)
  expect_equal(fitted(res.up), fitted(res.gen), tolerance = 1e-6)

  ## Each animal's breeding value, matched by its globulus id
  bv_up  <- ranef(res.up)$genetic
  bv_gen <- ranef(res.gen)$genetic
  orig   <- as.integer(names(bv_up)) - 100L
  expect_equal(as.numeric(bv_up)[match(as.integer(names(bv_gen)), orig)],
               as.numeric(bv_gen), tolerance = 1e-6)
})

test_that("the QC report of a recoded pedigree names animals by their ids (#58)", {

  ## preGSf90 reports animals by the codes of the pedigree breedR wrote, which
  ## are recoded here. The report used to show those codes, and the code
  ## of an excluded animal was another animal's id.
  N   <- max(globulus[, c('self', 'dad', 'mum')])
  flip <- function(x) ifelse(x == 0, 0L, as.integer(N + 1L - x))
  glob_rev <- globulus
  glob_rev[, c('self', 'dad', 'mum')] <-
    lapply(globulus[, c('self', 'dad', 'mum')], flip)

  ## Genotype some founders too, so that a genotyped animal has a genotyped
  ## parent. A local copy of the genotypes: the animal in row 3 fails the call
  ## rate, and the animal in row `kid` conflicts with its genotyped sire.
  founders <- setdiff(unique(c(globulus$dad, globulus$mum)),
                      c(0, globulus$self))
  ids <- c(globulus$self[1:50], founders[1:10])
  G   <- Gmat
  G[3, 1:150] <- 5L
  kid <- which(globulus$dad[1:50] %in% founders[1:10])[1]
  sire <- match(globulus$dad[kid], ids)
  G[kid, 151:200] <- 2L
  G[sire, 151:200] <- 0L

  qc_dir  <- breedR_workdir('qc_recoded_')
  on.exit(unlink(qc_dir, recursive = TRUE), add = TRUE)
  qc_file <- file.path(qc_dir, "test_ssgblup_qc_recoded.txt")
  write_snp_file(G, ids = flip(ids), file = qc_file)

  res.qc <- suppressWarnings(suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = glob_rev[, 1:3], id = 'self'),
            genomic = list(snp_file = qc_file, verify_parentage = 1L,
                           extra_options = 'thrStopCorAG -1'),
            data = glob_rev)))
  expect_false(is.null(attr(get_pedigree(res.qc), 'map')))

  ## the call-rate report names the animal of row 3 by its id
  qc <- res.qc$genomic
  expect_identical(qc$n_animals_excluded, 1L)
  expect_identical(qc$excluded_animals[[ncol(qc$excluded_animals) - 1L]],
                   flip(ids[3]))

  ## and the conflict names the progeny and its sire by their ids
  sire_line <- grep("Animal - Sire Conflict", qc$conflicts, value = TRUE)
  expect_length(sire_line, 1L)
  expect_identical(strsplit(trimws(sire_line), "\\s+")[[1]][5:6],
                   as.character(flip(ids[c(kid, sire)])))
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

test_that("gibbsf90(genomic=) builds the XrefID from a recoded pedigree (#55)", {

  ## Reversed codes put offspring before their parents, so the pedigree is
  ## recoded. gibbsf90() used to stop because it never passed the pedigree on.
  N   <- max(globulus[, c('self', 'dad', 'mum')])
  flip <- function(x) ifelse(x == 0, 0L, as.integer(N + 1L - x))
  glob_rev <- globulus
  glob_rev[, c('self', 'dad', 'mum')] <-
    lapply(globulus[, c('self', 'dad', 'mum')], flip)
  map <- attr(suppressWarnings(build_pedigree(1:3, data = glob_rev[, 1:3])),
              'map')
  expect_false(is.null(map))

  gdir <- breedR_workdir('gibbs_genomic_')
  on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
  f <- file.path(gdir, "test_gibbs_genomic.txt")
  write_snp_file(Gmat, ids = flip(gen_ids), file = f)

  gb <- suppressWarnings(suppressMessages(
    gibbsf90(phe_X ~ gg,
             genetic = list(model = 'add_animal',
                            pedigree = glob_rev[, 1:3], id = 'self'),
             genomic = list(snp_file = f, verify_parentage = 0L),
             data = glob_rev, n_samples = 100L, burnin = 10L, thin = 1L)))
  on.exit(unlink(gb$dir, recursive = TRUE), add = TRUE)

  ## each genotyped animal gets the code breedR wrote for its own id
  xref <- utils::read.table(file.path(gb$dir, paste0(basename(f), "_XrefID")))
  expect_equal(xref[[2]], flip(gen_ids))
  expect_equal(xref[[1]], map[flip(gen_ids)])

  expect_equal(gb$genomic$n_snp_total, nsnp)
  expect_true(any(grepl("Number of Genotyped Animals: 60", gb$output)))
  expect_gt(nrow(gb$solutions), 0)
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

test_that("postgsf90 works on a fit with OPTION missing (#121)", {
  ## A response spanning 0 makes progsf90() code missing records with
  ## OPTION missing, which postGSf90 rejects. Centring the response is
  ## absorbed by the fixed effects, so the SNP effects must not change.
  gl <- globulus
  gl$phe_c <- gl$phe_X - mean(gl$phe_X)
  fit <- function(resp, tag) {
    f <- file.path(tempdir(), paste0("test_gwas_missing_", tag, ".txt"))
    write_snp_file(Gmat, ids = gen_ids, file = f)
    file.remove(list.files(dirname(f), pattern = "_XrefID$",
                           full.names = TRUE))
    suppressMessages(
      remlf90(fixed = as.formula(paste(resp, "~ gg")),
              genetic = list(model = 'add_animal',
                             pedigree = gl[, 1:3], id = 'self'),
              genomic = list(snp_file = f, verify_parentage = 0L,
                             save_ginverse = TRUE),
              data = gl))
  }
  res.c <- fit("phe_c", "c")
  res.x <- fit("phe_X", "x")
  expect_true(any(grepl("^OPTION missing",
                        readLines(file.path(res.c$reml$dir, "parameters")))))

  ## snp_var also runs the BLUP pass, which keeps the option
  gwas.c <- postgsf90(res.c, snp_var = TRUE)
  gwas.x <- postgsf90(res.x)
  expect_equal(nrow(gwas.c$snp_sol), nsnp)
  expect_equal(gwas.c$snp_sol$solution, gwas.x$snp_sol$solution,
               tolerance = 1e-5)
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

test_that("predf90 computes reliabilities from postgsf90(snp_var = TRUE) (#119)", {

  rel_snp <- file.path(tempdir(), "test_predf90_acc_snp.txt")
  write_snp_file(Gmat, ids = gen_ids, file = rel_snp)
  file.remove(list.files(dirname(rel_snp), pattern = "_XrefID$",
                         full.names = TRUE))
  res.acc <- suppressMessages(
    remlf90(fixed = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = globulus[, 1:3], id = 'self'),
            genomic = list(snp_file = rel_snp, verify_parentage = 0L,
                           save_ginverse = TRUE),
            data = globulus))
  d <- res.acc$reml$dir
  sol_md5 <- tools::md5sum(file.path(d, "solutions"))

  gwas <- postgsf90(res.acc, snp_var = TRUE)

  ## The variances postGSf90 computed the SNP prediction error variances at
  ## are the REML estimates, not the starting values the fit's parameter file
  ## holds. With the starting values the reliabilities came out about three
  ## times too high, and nothing failed.
  par <- readLines(file.path(d, "parameters_postgs"))
  num <- function(head) as.numeric(par[which(par == head) + 1L])
  start <- as.numeric(readLines(file.path(d, "parameters"))[
    which(readLines(file.path(d, "parameters")) == "(CO)VARIANCES") + 1L])
  expect_equal(num("(CO)VARIANCES"), res.acc$var["genetic", 1],
               tolerance = 1e-6)
  expect_equal(num("RANDOM_RESIDUAL VALUES"), res.acc$var["Residual", 1],
               tolerance = 1e-6)
  expect_false(isTRUE(all.equal(start, res.acc$var["genetic", 1])))

  ## and so are those the inverse itself was computed at
  blup <- readLines(file.path(d, "parameters_blup"))
  expect_equal(as.numeric(blup[which(blup == "(CO)VARIANCES") + 1L]),
               res.acc$var["genetic", 1], tolerance = 1e-6)
  expect_equal(as.numeric(blup[which(blup == "RANDOM_RESIDUAL VALUES") + 1L]),
               res.acc$var["Residual", 1], tolerance = 1e-6)

  ## the fit's own solutions are left as they were
  expect_identical(tools::md5sum(file.path(d, "solutions")), sol_md5)
  ## and so are the SNP effects
  expect_equal(gwas$snp_sol$solution, postgsf90(res.acc)$snp_sol$solution)

  ## postgsf90(res.acc) above cleared them; ask again. One file per trait and
  ## correlated effect: snp_var_1_1 for this single-trait model.
  gwas <- postgsf90(res.acc, snp_var = TRUE)
  expect_identical(list.files(d, "^snp_var_"), "snp_var_1_1")

  new_snp <- file.path(tempdir(), "test_predf90_acc_new.txt")
  set.seed(11)
  write_snp_file(matrix(sample(0:2, 20 * nsnp, replace = TRUE,
                               prob = c(0.25, 0.5, 0.25)), 20),
                 ids = 900001:900020, file = new_snp)
  pred <- predf90(snp_file = new_snp, dir = gwas$dir, acc = TRUE)
  expect_equal(nrow(pred), 20L)
  expect_true("reliability" %in% names(pred))
  expect_true(all(is.finite(pred$reliability)))
  expect_true(all(pred$reliability > 0 & pred$reliability < 1))

  ## p-values, which needed the same inverse and so never worked
  pv <- postgsf90(res.acc, snp_p_value = TRUE)
  expect_equal(nrow(pv$pvalues), nsnp)
  expect_true(all(is.finite(pv$pvalues$neg_log10_pval)))

  ## A later call without them does not return the earlier ones
  plain <- postgsf90(res.acc)
  expect_null(plain$pvalues)
  expect_length(list.files(d, "^snp_var_"), 0L)

  ## and the option given the raw way works as well
  postgsf90(res.acc, extra_options = "snp_var")
  expect_length(list.files(d, "^snp_var_"), 1L)
})

test_that("predf90 errors without a prior postgsf90 run", {
  lone_snp <- file.path(tempdir(), "test_predf90_lone.txt")
  write_snp_file(Gmat, ids = gen_ids, file = lone_snp)
  expect_error(
    predf90(snp_file = lone_snp, dir = file.path(tempdir(), "no_postgs")),
    "snp_pred")
})
