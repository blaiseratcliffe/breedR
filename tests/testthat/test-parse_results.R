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


test_that("Parse a covariance block whose rows start flush with a minus", {

  ## Fortran carriage control usually leaves column 1 blank, but a negative
  ## value wide enough to fill its field starts flush against it. Reading the
  ## row must not consume that value as if it were the leading blank.
  flush_block <- c("  4.2000    -0.16000    ",
                   "-0.16000     0.50000    ")

  expect_error(m <- parse.txtmat(flush_block), NA)
  expect_equal(m, matrix(c(4.2, -0.16, -0.16, 0.5), 2, 2))
  expect_true(isSymmetric(m))

  ## Wider case: every row but the first starts flush.
  flush3 <- c(" 0.18966     -0.27594E-04 -0.60787E-04",
              "-0.27594E-04  0.39831E-02 -0.29742E-02",
              "-0.60787E-04 -0.29742E-02  0.29915E-02")
  expect_error(m3 <- parse.txtmat(flush3), NA)
  expect_identical(dim(m3), c(3L, 3L))
  expect_true(isSymmetric(m3))
  expect_equal(m3[2, 1], -0.27594e-04)

})


test_that("Parse the residual (co)variance echo of a real REML log", {

  ## The parameter echo near the top of the log is the one block the backend
  ## prints with negative values flush at column 1. remlf90(traits = ) reads it
  ## to work out which covariance parameters the AI matrix left free.
  x3 <- readLines(file.path(testdata, "airemlf90_log_3.txt"))
  at <- grep("Residual (co)variance Matrix", x3, fixed = TRUE)
  expect_identical(length(at), 1L)

  expect_error(m <- parse.txtmat(extract_block(at + 1L, x3)), NA)
  expect_equal(m, matrix(c(4.2, -0.16, -0.16, 0.5), 2, 2))

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
parse_hetres_fixture <- function(key, fixed, random = NULL, var.ini = list(),
                                 group_name = 'g') {
  data <- data.frame(x  = seq(0, 2, length.out = 100),
                     g  = factor(rep(1:50, 2)),
                     y  = rep(c(9, 11), 50),
                     y3 = rep(c(2, 4), 50))
  names(data)[names(data) == 'g'] <- group_name
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


test_that("hetres standard errors are unaffected by random-effect names", {

  control <- parse_hetres_fixture(1, y ~ x, ~ g, list(g = 3.4))
  for (nm in c('a0', 'a1')) {
    ## Renaming the group changes no estimates, but duplicates a coefficient
    ## name in invAI. A name lookup would select the group's S.E. instead.
    res <- parse_hetres_fixture(1, y ~ x, reformulate(nm),
                                setNames(list(3.4), nm), group_name = nm)
    expect_identical(rownames(res$var), nm)
    expect_equal(res$hetres, control$hetres)
  }
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


test_that("parse_results() surfaces the backend's own diagnostic when the run never starts (issue #14)", {

  ## BLUPF90+ can exit 0 after a fatal, setup-time error -- a missing data
  ## file, a malformed NUMBER_OF_EFFECTS, a truncated parameter file -- see
  ## #14. None of those reach the REML loop, so unlike a run that merely
  ## fails to converge (tested above), the log has no 'In round' line at all.
  ## Before the fix, parse_results() read past this straight into
  ## last.round.idx <- tail(grep('In round', reml.out), 1), which comes back
  ## integer(0), and died several lines later on reml.out[last.round.idx]
  ## with an opaque "subscript out of bounds" instead of the diagnostic that
  ## was sitting in reml.out the whole time.
  data <- data.frame(x = seq(0, 2, length.out = 100),
                     g = factor(rep(1:50, 2)),
                     y = rep(c(9, 11), 50))
  mc <- call('remlf90', fixed = quote(y ~ x), random = quote(~ g), data = quote(data))
  mf <- build.mf(mc)
  effects <- build.effects(mf, NULL, NULL, NULL, list(g = 3.4))

  real_log <- readLines(file.path(testdata, 'airemlf90_log_hetres_1.txt'))
  no_rounds <- c(real_log[!grepl('In round', real_log)],
                "There is no such data file: no_such_file")

  expect_error(
    parse_results(file.path(testdata, 'airemlf90_sol_hetres_1.txt'),
                  effects, mf, no_rounds, 'ai', quote(remlf90())),
    "no such data file", fixed = TRUE
  )
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


# Captured BLUPF90+ 2.76 analytical fit: the omitted random effect is y2,
# while the retained y1 has variance 79/15 and BLUPs (79/84)*c(-3,-2,-1,1,2,3).
parse_issue51_fixture <- function(method = 'ai', nested = FALSE,
                                  sol = NULL, out = NULL) {
  data <- data.frame(rep = factor(rep(1:6, each = 4)),
                     y1 = 10 + rep(c(-3,-2,-1,1,2,3), each = 4) +
                       rep(c(-1,1,-1,1), 6),
                     y2 = 20 + rep(c(-2,-2,2,2), 6))
  mc <- call('remlf90', fixed = cbind(y1, y2) ~ 1,
             random = if (!nested) ~ rep, data = quote(data))
  mf <- build.mf(mc)
  if (nested) {
    incidence <- matrix(0, 24, 6)
    g <- as.integer(data$rep)
    incidence[cbind(1:24, g)] <- .7
    incidence[cbind(1:24, g %% 6 + 1)] <- .3
    generic <- list(rep = list(incidence = incidence, covariance = diag(6),
                               var.ini = diag(c(5, 0))))
  } else generic <- NULL
  effects <- build.effects(mf, NULL, NULL, generic,
                            list(rep = diag(c(5, 0))),
                            traits = list(rep = c(y1 = TRUE, y2 = FALSE)))
  key <- paste0('issue51_', if (nested) 'nested' else 'analytic', '_', method)
  if (is.null(sol)) sol <- file.path(testdata, paste0(key, '.sol'))
  if (is.null(out)) out <- readLines(file.path(testdata, paste0(key, '.log')))
  parse_results(sol, effects, mf, out, method, quote(remlf90()))
}


test_that('restricted AI and EM results keep full trait coordinates', {
  for (method in c('ai', 'em')) {
    ans <- parse_issue51_fixture(method)
    expect_identical(names(ans$ranef$rep), c('y1', 'y2'))
    expect_equal(ans$ranef$rep$y1$value, (79/84) * c(-3,-2,-1,1,2,3),
                 tolerance = 1e-5)
    expect_true(all(is.na(as.matrix(ans$ranef$rep$y2))))
    expect_identical(rownames(ans$ranef$rep$y1), as.character(1:6))
    cv <- if (method == 'ai') ans$var[, 1] else ans$var
    expect_equal(cv$rep['y1', 'y1'], 5.2667)
    expect_true(all(is.na(cv$rep[2, ])))
    expect_true(all(is.na(cv$rep[, 2])))
    expect_identical(dim(cv$rep), c(2L, 2L))
    expect_equal(unname(diag(cv$Residual)), c(1.3333, 4.1739))
    expect_identical(cv$Residual[1, 2], 0)
    if (method == 'ai') {
      expect_true(all(is.na(ans$var[['rep', 'S.E.']][2, ])))
      expect_identical(rownames(ans$reml$invAI),
                       c('rep.y1', 'resid.y1', 'resid.y2'))
    } else expect_null(ans$reml$invAI)
  }
})


test_that('restricted solutions align explicit trait, effect and level keys', {
  path <- file.path(testdata, 'issue51_analytic_ai.sol')
  rows <- read.table(path, skip = 1)
  write_sol <- function(x) {
    f <- tempfile()
    writeLines('trait effect level solution se', f)
    write.table(x, f, append = TRUE, row.names = FALSE, col.names = FALSE)
    f
  }
  reference <- parse_issue51_fixture()
  shuffled <- parse_issue51_fixture(sol = write_sol(rows[nrow(rows):1, ]))
  expect_identical(shuffled$ranef, reference$ranef)
  expect_identical(shuffled$fixed, reference$fixed)
  retained <- rows[!(rows$V1 == 2 & rows$V2 == 2), ]
  sparse <- parse_issue51_fixture(sol = write_sol(retained))
  expect_identical(sparse$ranef, reference$ranef)
  expect_error(parse_issue51_fixture(sol = write_sol(rbind(rows, rows[1, ]))),
               'duplicate effect/trait/level')
  expect_error(parse_issue51_fixture(sol = write_sol(rows[-1, ])),
               'missing active solution')
  bad <- rows
  bad$V3[1] <- 99
  expect_error(parse_issue51_fixture(sol = write_sol(bad)), 'out-of-range')
  bad <- rows
  bad$V1[1] <- 3
  expect_error(parse_issue51_fixture(sol = write_sol(bad)), 'out-of-range')
})


test_that('restricted nested effects use their positive-level solution anchor', {
  ans <- parse_issue51_fixture(nested = TRUE)
  lay <- pf90_effect_layout(ans$effects, 2)
  expect_identical(as.integer(lay$effect), c(1L, 3L))
  expect_identical(names(ans$ranef), 'rep')
  expect_identical(dim(ans$ranef$rep$y1), c(6L, 2L))
  expect_true(all(is.finite(ans$ranef$rep$y1$value)))
  expect_true(all(is.na(ans$ranef$rep$y2$value)))
})


test_that('restricted nonconvergence preserves covariance dimensions and names', {
  out <- readLines(file.path(testdata, 'issue51_analytic_ai.log'))
  at <- tail(grep('In round', out), 1)
  out[at] <- sub('In round\\s+[0-9]+',
                 paste('In round', MAX_REML_ITERATIONS), out[at])
  expect_warning(ans <- parse_issue51_fixture(out = out), 'did not converge')
  expect_identical(dim(ans$var), c(2L, 2L))
  expect_identical(dim(ans$var[['rep', 1]]), c(2L, 2L))
  expect_identical(rownames(ans$var[['rep', 1]]), c('y1', 'y2'))
  expect_true(all(is.na(ans$var[['rep', 1]])))
  expect_identical(rownames(ans$reml$invAI),
                   c('rep.y1', 'resid.y1', 'resid.y2'))
})


test_that('restricted AI dimension mismatches fail explicitly', {
  out <- readLines(file.path(testdata, 'issue51_analytic_ai.log'))
  at <- grep('Residual (co)variance Matrix', out, fixed = TRUE)
  out[at + 1:2] <- c(' 1 .2', ' .2 4')
  expect_error(parse_issue51_fixture(out = out),
               'Backend AI layout does not match the active covariance parameters')
})


test_that('touching fixed-width solution keys retain original trait numbers', {
  f <- tempfile()
  writeLines(c('trait/effect level solution se',
                '   1   1         1        10.0        1.0',
                '   2   1         1        20.0        1.0',
                '   1   1         2        11.0        1.0',
                '   2   1         2        21.0        1.0',
                '   11000         1         3.0        1.0',
                '   21000         1         0.0        0.0'), f)
  ans <- read_pf90_solutions(f, restricted = TRUE)
  expect_equal(ans$trait, rep(1:2, 3))
  expect_equal(ans$effect, c(1,1,1,1,1000,1000))
  expect_equal(ans$value, c(10,20,11,21,3,0))
  # A model may contain only a nested random effect, so all rows can merge.
  writeLines(c('trait/effect level solution se',
                '   11000         1         3.0        1.0',
                '   21000         1         0.0        0.0'), f)
  ans <- read_pf90_solutions(f, restricted = TRUE)
  expect_equal(ans$trait, 1:2)
  expect_equal(ans$effect, c(1000, 1000))
})
