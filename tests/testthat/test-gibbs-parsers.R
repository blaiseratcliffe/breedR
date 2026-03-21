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

test_that("parse_gibbs_results reads gibbs_samples", {
  dir <- file.path(tempdir(), "test_gibbs_parse2")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))

  # Mock gibbs_samples: 3 header lines, then alternating round/values
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
