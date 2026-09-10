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


test_that("first_number() treats a missing line as no number", {

  ## Leading number, hence the name: the caller splits the log line first and
  ## feeds the remainders, because "-2logL = 39.76" leads with the -2.
  expect_equal(first_number("  39.76  logL convergence 0.1E-06"), 39.76)
  expect_equal(first_number("-2logL =  39.76 : AIC = 43.76"), -2)
  expect_equal(first_number(c("x 1", "no digits here")), c(1, NA))

  ## NA in, NA out. regexpr() answers NA for an NA input, and using that as a
  ## subscript used to fail with 'replacement has length zero' -- which is how
  ## a log line that simply was not there got reported.
  expect_identical(first_number(NA_character_), NA_real_)
  expect_equal(first_number(c("v 2", NA)), c(2, NA))
  expect_identical(first_number(character(0)), numeric(0))
})


test_that("the log-likelihood is found when the backend bends the AI matrix", {

  ## BLUPF90+ puts a 'Corrections made ... bending proportions' line between
  ## the -2logL line and the round it belongs to whenever it has to bend the
  ## AI matrix -- routine for a negative rho. parse_results() used to take the
  ## -2logL line by a fixed offset from the last round, so a bent fit read the
  ## bending line, got NA, and died several frames later. It is searched for
  ## now, so both shapes give the same answer.
  plain <- c("-2logL =     39.7649 : AIC =     43.7649  logL convergence 0.1E-06",
             "  In round          183  convergence=  2.9E-005",
             "  delta convergence=  9.9E-007")
  bent  <- append(plain,
                  "Corrections made:    1 , final bending proportions of AI and EM",
                  after = 1)

  logl_of <- function(out) {
    last.round.idx <- tail(grep('In round', out), 1)
    idx <- tail(grep('-2logL', out[seq_len(last.round.idx)]), 1)
    first_number(strsplit(strsplit(out[idx], split = '-2logL =')[[1]][2],
                          split = ': AIC =')[[1]])
  }

  expect_equal(logl_of(plain)[1], 39.7649)
  expect_equal(logl_of(bent)[1],  39.7649)
})


## Parse a heterogeneous-residual fixture pair (see inst/testdata/index.json).
## parse_results() takes the effect structure, the factor levels and the trait
## names from mf and effects, never the data values, so a small deterministic
## data frame of the same shape stands in for the simulated data of the run.
parse_hetres_fixture <- function(key, fixed, random = NULL, var.ini = list()) {
  data <- data.frame(x  = seq(0, 2, length.out = 100),
                     g  = factor(rep(1:50, 2)),
                     y  = rep(c(9, 11), 50),
                     y3 = rep(c(2, 4), 50))
  mc <- call('remlf90', fixed = fixed, random = random, data = quote(data))
  mf <- build.mf(mc)
  effects <- build.effects(mf, NULL, NULL, NULL, var.ini)
  parse_results(
    file.path(testdata, paste0('airemlf90_sol_hetres_', key, '.txt')),
    effects, mf,
    readLines(file.path(testdata, paste0('airemlf90_log_hetres_', key, '.txt'))),
    'ai', quote(remlf90())
  )
}


test_that("parse_results() reads a heterogeneous residual variance fit (issue #1)", {

  ## One trait, y ~ x + random g (50 levels), log(var(e)) = a0 + a1*x.
  ## Under hetres the backend prints the coefficients in place of a residual
  ## variance block, which used to trip the count of variance components.
  expect_error(
    res <- parse_hetres_fixture(1, y ~ x, ~ g, list(g = 3.4)),
    NA
  )

  ## no Residual row: the coefficients replace it
  expect_identical(rownames(res$var), 'g')
  expect_equal(res$var['g', 'Estimated variances'], 2.0684)
  expect_equal(res$var['g', 'S.E.'], 0.43550)

  expect_identical(rownames(res$hetres), c('a0', 'a1'))
  expect_identical(colnames(res$hetres), c('Estimate', 'S.E.'))
  expect_equal(unname(res$hetres[, 'Estimate']),
               c(0.576236467472488, 0.779007867141735))

  ## S.E. from the inverse AI matrix, which follows the G components.
  ## The backend's own 'SE for R' gives a0's only, and agrees with it.
  expect_equal(unname(res$hetres[, 'S.E.']), sqrt(c(0.39831E-02, 0.29915E-02)))
  expect_equal(unname(res$hetres['a0', 'S.E.']), 0.063112, tolerance = 1e-4)
  expect_identical(dimnames(res$reml$invAI)[[1]], c('g', 'a0', 'a1'))

  expect_equal(res$fit$'-2logL', 8530.66078354080)
  expect_length(res$funvars, 0L)

  ## the location effects are still read from the solutions file
  expect_equal(res$fixed$x[[1]]$value, 1.07301691)
  expect_identical(rownames(res$ranef$g[[1]]), as.character(1:50))
})


test_that("hetres coefficients of a multi-trait fit keep their trait order", {

  ## Two traits, no random effect, one covariate per trait. The backend prints
  ## the coefficients coefficient-major, trait-inner, and the inverse AI matrix
  ## follows the same order. A trait-major mapping would give a0.y3 the S.E.
  ## of a1.y.
  expect_error(
    res <- parse_hetres_fixture(2, cbind(y, y3) ~ x),
    NA
  )

  expect_identical(rownames(res$hetres), c('a0.y', 'a0.y3', 'a1.y', 'a1.y3'))
  expect_equal(unname(res$hetres[, 'Estimate']),
               c(0.448080563362213, 2.02177454365617,
                 0.870952863693744, -1.19361378589017))
  expect_equal(unname(res$hetres[, 'S.E.']),
               sqrt(c(0.26417E-02, 0.26226E-02, 0.19432E-02, 0.19251E-02)))
  expect_identical(dimnames(res$reml$invAI)[[1]], rownames(res$hetres))

  expect_false('Residual' %in% rownames(res$var))
  expect_identical(nrow(res$var), 0L)
})


test_that("homoscedastic logs are not taken for hetres output", {

  ## Detection relies on lines that only hetres output contains, so the
  ## parsing of every other log is unchanged.
  for (f in c(paste0('airemlf90_log_', 1:4, '.txt'), 'remlf90_log_1.txt'))
    expect_false(any(grepl(hetres_coef_re, readLines(file.path(testdata, f)))),
                 label = f)
})


test_that("parse_hetres() refuses coefficients in an unexpected order", {

  ok <- c(" new R",
          "           1 -th trait:           1 -th coefficient =  0.4",
          "           2 -th trait:           1 -th coefficient =  2.0",
          "           1 -th trait:           2 -th coefficient =  0.8",
          "           2 -th trait:           2 -th coefficient = -1.2",
          " inverse of AI matrix (Sampling Variance)")
  expect_equal(unname(parse_hetres(ok, 2, c('y', 'y3'))[, 'Estimate']),
               c(0.4, 2.0, 0.8, -1.2))

  ## trait-major printing, or a count that does not fill every trait, must
  ## fail rather than be assigned to the wrong coefficients
  swapped <- ok[c(1, 2, 4, 3, 5, 6)]
  expect_error(parse_hetres(swapped, 2, c('y', 'y3')), 'hetres')
  expect_error(parse_hetres(ok[-5], 2, c('y', 'y3')), 'hetres')
})
