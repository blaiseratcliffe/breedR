context("parse_results")

test_that("Identify blocks of large covariance matrices", {
  
  ## airemlf90_log_4.txt - 10-variate model.
  ## with large residual covariance matrix (10x10) which wraps lines
  reml_log <- readLines(file.path(testdata, "airemlf90_log_4.txt"))
  resid_var_ini_line <- grep("Residual variance", reml_log)
  
  expect_error(
    bl <- extract_block(resid_var_ini_line + 1, reml_log),
    NA
  )
  
  expect_identical(length(bl), 20L)
  
})

test_that("Parse large covariance matrices", {
  
  test_line <- c("   1.23  4.56   7.89E-02",
                " 1 2 3",
                " 4")
  
  non_square_test <- rep(test_line, 2)
  
  expect_error(
    parse.txtmat(non_square_test),
    "square matrix"
  )
  expect_error(
    ns_mat <- parse.txtmat(non_square_test, square = FALSE),
    NA
  )
  expect_identical(dim(ns_mat), c(2L, 7L))

  
  square_test <- rep(test_line, 7)
  
  expect_error(
    s_mat <- parse.txtmat(square_test),
    NA
  )
  expect_identical(dim(s_mat), c(7L, 7L))


})


test_that("Parse the fit line when fields carry trailing text", {

  ## current BLUPF90+ appends the logL convergence to the AIC field
  fit_line <- paste("-2logL =     5748.87709756475     : AIC =",
                    "    5754.87709756475    logL convergence",
                    "  0.727737843402541E-14")

  fields <- strsplit(strsplit(fit_line, split = '-2logL =')[[1]][2],
                     split = ': AIC =')[[1]]

  ## coercing the whole field loses the AIC
  expect_true(is.na(suppressWarnings(as.numeric(fields[2]))))

  expect_warning(ans <- first_number(fields), NA)
  expect_equal(ans, c(5748.87709756475, 5754.87709756475))

  ## older output, with nothing trailing the AIC
  old_line <- "-2logL =     5748.87709756475     : AIC =     5754.87709756475"
  expect_equal(
    first_number(strsplit(strsplit(old_line, split = '-2logL =')[[1]][2],
                          split = ': AIC =')[[1]]),
    c(5748.87709756475, 5754.87709756475)
  )

  ## exponents, signs and absent numbers
  expect_equal(first_number("  2.546803377662850E-014"), 2.54680337766285e-14)
  expect_equal(first_number(" -1234.5 then words"), -1234.5)
  expect_equal(first_number("42"), 42)
  expect_identical(first_number("no number here"), NA_real_)
  expect_identical(first_number(character(0)), numeric(0))
})
