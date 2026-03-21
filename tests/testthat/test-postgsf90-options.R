## Tests for PostGSF90 option building (R/genomic.R)

context("PostGSF90 option building")

test_that("windows_variance generates correct option", {
  opts <- build_postgsf90_options(
    windows_variance = 20, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("windows_variance 20" %in% opts)
  expect_true("windows_variance_type 1" %in% opts)
})

test_that("windows_variance_mbp generates correct option", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = 1,
    windows_variance_type = 2, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("windows_variance_mbp 1" %in% opts)
  expect_true("windows_variance_type 2" %in% opts)
})

test_that("no windows options when both NULL", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_false(any(grepl("windows_variance", opts)))
})

test_that("manhattan_plot generates option", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = TRUE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("Manhattan_plot_R" %in% opts)
})

test_that("which_weight generates option", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = 1,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("which_weight 1" %in% opts)
})

test_that("snp_p_value generates option", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = TRUE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("snp_p_value" %in% opts)
})

test_that("postgs_trt_eff generates option with two values", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = c(1, 2), extra_options = NULL)
  expect_true("postgs_trt_eff 1 2" %in% opts)
})

test_that("postgs_trt_eff errors on wrong length", {
  expect_error(
    build_postgsf90_options(
      windows_variance = NULL, windows_variance_mbp = NULL,
      windows_variance_type = 1, which_weight = NULL,
      manhattan_plot = FALSE, snp_p_value = FALSE,
      postgs_trt_eff = c(1), extra_options = NULL))
})

test_that("extra_options passed through", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = c("solutions_postGS mysol"))
  expect_true("solutions_postGS mysol" %in% opts)
})

test_that("empty options returns empty vector", {
  opts <- build_postgsf90_options(
    windows_variance = NULL, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = NULL,
    manhattan_plot = FALSE, snp_p_value = FALSE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_length(opts, 0)
})

test_that("combined options all present", {
  opts <- build_postgsf90_options(
    windows_variance = 20, windows_variance_mbp = NULL,
    windows_variance_type = 1, which_weight = 2,
    manhattan_plot = TRUE, snp_p_value = TRUE,
    postgs_trt_eff = NULL, extra_options = NULL)
  expect_true("windows_variance 20" %in% opts)
  expect_true("which_weight 2" %in% opts)
  expect_true("Manhattan_plot_R" %in% opts)
  expect_true("snp_p_value" %in% opts)
})
