## Tests for Gibbs result parsers (R/gibbs.R)
## Uses mock files — no GIBBSF90+/POSTGIBBSF90 binary needed.

context("Gibbs result and diagnostics parsers")

## -- parse_gibbs_results() --

test_that("parse_gibbs_results reads final_solutions", {
  dir <- file.path(tempdir(), "test_gibbs_parse1")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))

  writeLines(c(
    "trait/effect level  solution        SD",
    "   1   1         1         13.58          0.50",
    "   1   1         2         14.12          0.82",
    "   1   2         1          2.72          2.30"
  ), file.path(dir, "final_solutions"))

  result <- parse_gibbs_results(dir)
  expect_equal(nrow(result$solutions), 3)
  expect_equal(names(result$solutions), c("trait", "effect", "level",
                                            "solution", "sd"))
  expect_equal(result$solutions$solution[1], 13.58)
})

## GIBBSF90+ opens gibbs_samples with one line giving the number of
## (co)variance components, then one mapping line per component, then
## alternating round and value lines. The header is therefore as long as the
## model is wide, and only a model with exactly two components has the three
## header lines the parser used to assume (#77). The samples below are real
## output: a fixed-effects model, one and two random effects, and two traits
## with one random effect (3 genetic + 3 residual components).
samples_dir <- function(lines, key) {
  dir <- file.path(tempdir(), paste0("test_gibbs_samples_", key))
  dir.create(dir, showWarnings = FALSE)
  writeLines(lines, file.path(dir, "gibbs_samples"))
  dir
}

test_that("parse_gibbs_results reads samples of a one-component model", {
  dir <- samples_dir(c(
    "      -1       1       4",
    "    1    0    0    1    1",
    "      10       1",
    "   10.27    ",
    "      20       1",
    "   9.967    "
  ), "1c")
  on.exit(unlink(dir, recursive = TRUE))

  samples <- parse_gibbs_results(dir)$samples
  expect_equal(dim(samples), c(2L, 1L))
  expect_equal(unname(samples[, 1]), c(10.27, 9.967))
  expect_equal(rownames(samples), c("10", "20"))
})

test_that("parse_gibbs_results reads samples of a three-component model", {
  dir <- samples_dir(c(
    "      -1       3       4",
    "    1    3    3    1    1",
    "    2    4    4    1    1",
    "    3    0    0    1    1",
    "      10       3",
    "  0.5074E-01   2.123       8.635    ",
    "      20       3",
    "  0.8970E-01   3.680       8.668    "
  ), "3c")
  on.exit(unlink(dir, recursive = TRUE))

  samples <- parse_gibbs_results(dir)$samples
  expect_equal(dim(samples), c(2L, 3L))
  expect_equal(unname(samples[1, ]), c(0.5074e-01, 2.123, 8.635))
  expect_equal(unname(samples[2, 3]), 8.668)
  expect_equal(rownames(samples), c("10", "20"))
})

test_that("parse_gibbs_results reads samples of a multi-trait model", {
  dir <- samples_dir(c(
    "      -1       6       4",
    "    1    3    3    1    1",
    "    2    3    3    1    2",
    "    3    3    3    2    2",
    "    4    0    0    1    1",
    "    5    0    0    1    2",
    "    6    0    0    2    2",
    "      10       6",
    paste("  0.1013     -0.2511       2.442       10.61",
          "     0.2098       3.751    ")
  ), "6c")
  on.exit(unlink(dir, recursive = TRUE))

  samples <- parse_gibbs_results(dir)$samples
  expect_equal(dim(samples), c(1L, 6L))
  expect_equal(unname(samples[1, ]),
               c(0.1013, -0.2511, 2.442, 10.61, 0.2098, 3.751))
})

test_that("parse_gibbs_results reads samples that start at a round line", {
  ## The header is found rather than assumed, so a file without one parses
  ## instead of collapsing to nothing.
  dir <- samples_dir(c(
    "      10       1",
    "   10.27    ",
    "      20       1",
    "   9.967    "
  ), "nohdr")
  on.exit(unlink(dir, recursive = TRUE))

  samples <- parse_gibbs_results(dir)$samples
  expect_equal(dim(samples), c(2L, 1L))
  expect_equal(unname(samples[, 1]), c(10.27, 9.967))
})

test_that("parse_gibbs_results reads gibbs_samples", {
  dir <- file.path(tempdir(), "test_gibbs_parse2")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))

  # Mock gibbs_samples: two components, so 3 header lines, then alternating
  # round/values
  writeLines(c(
    "      -1       2       4",
    "    1    2    2    1    1",
    "    2    0    0    1    1",
    "     100       2",
    "   3.500       14.200",
    "     200       2",
    "   3.800       13.900",
    "     300       2",
    "   4.100       14.500"
  ), file.path(dir, "gibbs_samples"))

  result <- parse_gibbs_results(dir)
  expect_equal(nrow(result$samples), 3)
  expect_equal(ncol(result$samples), 2)
  expect_equal(result$samples[1, 1], 3.5)
  expect_equal(result$samples[3, 2], 14.5)
  expect_equal(rownames(result$samples), c("100", "200", "300"))
})

test_that("parse_gibbs_results reads fort.99", {
  dir <- file.path(tempdir(), "test_gibbs_parse3")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))

  writeLines(c("5800.5", "5790.2", "5785.1"),
             file.path(dir, "fort.99"))

  result <- parse_gibbs_results(dir)
  expect_length(result$deviance, 3)
  expect_equal(result$deviance[1], 5800.5)
})

test_that("parse_gibbs_results handles missing files gracefully", {
  dir <- file.path(tempdir(), "test_gibbs_parse4")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))

  result <- parse_gibbs_results(dir)
  expect_null(result$solutions)
  expect_null(result$samples)
  expect_null(result$deviance)
})

## -- the gibbsf90() result --

test_that("gibbsf90() returns the pedigree that translates its genetic levels (#64)", {

  ## Animal 1 precedes its parents 2 and 3, so the pedigree is recoded with
  ## map 3 1 2: each code is also another animal's id
  ped <- data.frame(self = 1:3, dad = c(2L, 0L, 0L), mum = c(3L, 0L, 0L))
  dat <- data.frame(ped, y = c(1.2, 0.3, 2.1))
  recoded <- suppressWarnings(build_pedigree(1:3, data = ped))
  expect_identical(attr(recoded, 'map'), c(3L, 1L, 2L))

  bin <- breedR_workdir()
  on.exit(unlink(bin, recursive = TRUE), add = TRUE)
  gibbs_name <- if (breedR.os.type() == 'windows') 'gibbsf90+.exe' else 'gibbsf90+'
  file.create(file.path(bin, gibbs_name))

  ## A stand-in for gibbsf90+. gibbsf90() runs it with base::system2(), which
  ## can only be mocked in base; every other command goes to the real one.
  ## Like the binary, it reports the genetic effect (effect 2) by the codes
  ## breedR wrote in the data file (column 3), giving each level the phenotype
  ## (column 1) of the animal it read under that code.
  real_system2 <- base::system2
  local_mocked_bindings(
    system2 = function(command, ...) {
      if (basename(command) != gibbs_name) return(real_system2(command, ...))
      d <- utils::read.table('data')
      writeLines(c('trait/effect level  solution        SD',
                   paste(1L, 2L, d[[3]], d[[1]], 0)),
                 'final_solutions')
      character(0)
    }, .package = 'base')

  res <- suppressWarnings(
    gibbsf90(y ~ 1,
             genetic = list(model = 'add_animal', pedigree = ped, id = 'self'),
             data = dat, breedR.bin = bin, n_samples = 10L))
  on.exit(unlink(res$dir, recursive = TRUE), add = TRUE)

  expect_identical(res$pedigree, recoded)

  ## through the map, each animal gets its own value back
  sol <- res$solutions[res$solutions$effect == 2, ]
  id  <- match(sol$level, attr(res$pedigree, 'map'))
  expect_equal(sol$solution[match(dat$self, id)], dat$y)
  ## with the level read as the id, each animal gets another animal's value
  expect_false(isTRUE(all.equal(sol$solution[match(dat$self, sol$level)],
                                dat$y)))
})

## -- parse_postout_tables() --

test_that("parse_postout_tables extracts MCE table", {
  lines <- c(
    "                ********      Monte Carlo Error by Time Series      ********",
    "Pos.  eff1    eff2    trt1    trt2         MCE      Mean        HPD        Effective  Median   Mode  Independent",
    "                                                           Interval (95%)  sample size               chain size",
    "1     2       2       1       1     0.232      3.931   0.602     7.207     59.3     3.759  3.788           92",
    "2     0       0       1       1     0.205     14.122  11.020    17.490     65.3    14.190 14.005           92",
    "                ********      Posterior Standard Deviation      ********",
    "Pos.  eff1    eff2    trt1    trt2        PSD       Mean       PSD         Convergence     Auto-correlations   Independent",
    "1     2       2       1       1     1.791      3.931   0.421     7.441   -0.13      0.920  0.487  0.059          19",
    "2     0       0       1       1     1.656     14.122  10.876    17.369    0.14      0.746  0.410  0.051          19"
  )

  result <- parse_postout_tables(lines)

  # MCE table
  expect_false(is.null(result$mce_table))
  expect_equal(nrow(result$mce_table), 2)
  expect_equal(result$mce_table$mean[1], 3.931)
  expect_equal(result$mce_table$mean[2], 14.122)

  # HPD
  expect_false(is.null(result$hpd))
  expect_equal(result$hpd$hpd_lower[1], 0.602)
  expect_equal(result$hpd$hpd_upper[2], 17.490)

  # Effective size
  expect_equal(result$effective_size, c(59.3, 65.3))

  # PSD table
  expect_false(is.null(result$psd_table))
  expect_equal(nrow(result$psd_table), 2)

  # Geweke
  expect_equal(result$geweke, c(-0.13, 0.14))

  # Autocorrelations
  expect_false(is.null(result$autocorrelations))
  expect_equal(result$autocorrelations$autocorr_lag1[1], 0.920)
  expect_equal(result$autocorrelations$autocorr_lag50[2], 0.051)
})

test_that("parse_postout_tables handles empty input", {
  result <- parse_postout_tables(character(0))
  expect_length(result, 0)
})

test_that("parse_postout_tables handles no data lines", {
  lines <- c("some header text", "no numeric data here")
  result <- parse_postout_tables(lines)
  expect_length(result, 0)
})
