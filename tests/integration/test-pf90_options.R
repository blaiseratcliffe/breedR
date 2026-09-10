suppressPackageStartupMessages(
  require(spam)
)

### Test the interface to PROGSF90 OPTIONS ###

context("PROGSF90 options")

test_that('remlf90 parses progsf90.options correctly', {

  N = 100
  dat <- transform(data.frame(x = runif(N)),
                    y = 1 + 2*x + rnorm(N))
  
  ## Returns the fit, so the parameter-file checks below can read the file
  ## this particular run wrote: each fit works in its own directory, and
  ## tempdir() no longer holds anybody's 'parameters'.
  expect_opt <- function(opt, regexp) {
    res <- suppressMessages(
      remlf90(y ~ x,
              data = dat,
              progsf90.options = opt)
    )

    for (i in seq_along(opt))
    expect_true(any(grepl(regexp[i], res$reml$output)),
                label = paste('option', opt[i], 'passed correctly'))

    invisible(res)
  }

  ## some additional option
  expect_opt('tol 1d-01', 'tolerance .*? 0\\.1')
  ## BLUPF90+ acknowledges this one as
  ##   * Run EM-REML for 2 rounds and switch to AI-REML
  ## The assertion used to look for 'EM-REML iterations 2', which the backend
  ## has never printed, so it had been failing since it was written.
  expect_opt('EM-REML 2', 'Run EM-REML for\\s+2 rounds')

  ## Conflicting option: sol se
  # included
  res <- expect_opt('sol se', 'store solutions and s\\.e\\.')
  # only once
  parameters_file <- readLines(file.path(res$reml$dir, 'parameters'))
  expect_identical(length(grep('OPTION sol se', parameters_file)), 1L)

  ## Conflicting option: missing
  # if explicit, the one set by the user is used
  res <- expect_opt('missing 12345', 'missing observation .*? 12345')
  # only once
  parameters_file <- readLines(file.path(res$reml$dir, 'parameters'))
  expect_identical(length(grep('OPTION missing', parameters_file)), 1L)

  ## Multiple options, some conflicting, some not
  opts <- c('sol se', 'missing 12345', 'tol 1d-01', 'EM-REML 2')
  expr <- c('store solutions and s\\.e\\.',
            'missing observation .*? 12345',
            'tolerance .*? 0\\.1',
            'Run EM-REML for\\s+2 rounds')
  res <- expect_opt(opts, expr)
  parameters_file <- readLines(file.path(res$reml$dir, 'parameters'))

  ## remlf90() adds 'method VCE' of its own, so the file always carries one
  ## OPTION more than the user asked for and a bare count against
  ## length(opts) could never hold -- which is why this had been failing.
  ## Assert what the test is actually for instead: each requested option is
  ## written, written once (the point of the conflict handling above), and
  ## nothing else is in there.
  option_lines <- grep('^OPTION ', parameters_file, value = TRUE)

  for (o in opts)
    expect_identical(sum(option_lines == paste('OPTION', o)), 1L,
                     label = paste('option', o, 'written exactly once'))

  expect_setequal(option_lines, paste('OPTION', c(opts, 'method VCE')))
})


test_that('AI-remlf90() returns heritability and inverse AI matrix', {
  
  ## Simulate a dataset with a heritability of 2/(1+2+1+1) = 0.4
  set.seed(1234)
  dat <- breedR.sample.phenotype(fixed   = c(mu = 10, x = 2),
                                 random = list(u = list(nlevels = 3,
                                                        sigma2  = 1)),
                                 genetic = list(model    = 'add_animal',
                                                Nparents = c(10, 10),
                                                sigma2_a = 2,
                                                check.factorial = FALSE),
                                 spatial = list(model     = 'AR',
                                                grid.size = c(15, 15),
                                                rho       = c(.7, .8),
                                                sigma2_s  = 1),
                                 residual.variance = 1)
  
  res <- remlf90(phenotype ~ 1 + X.x,
                 random = ~ u,
                 genetic = list(model = 'add_animal',
                                pedigree = dat[, 1:3],
                                id = 'self'),
                 spatial = list(model = 'AR',
                                coord = dat[, c('Var1', 'Var2')],
                                rho   = c(.7, .8)),
                 data   = dat)

  # AIREMLF90 output
  expect_true(any(grepl('* SE for function of \\(co\\)variances Heritability', 
                        res$reml$output)))
  expect_true(any(grepl('Heritability  - Function: ', res$reml$output)))
  
  # parsed heritability and inverse AI matrix
  expect_is(res$funvars, 'matrix')
  expect_identical(rownames(res$funvars),
                   c('mean', 'sample mean', 'sample sd'))
  expect_is(res$reml$invAI, 'matrix')
  expect_identical(dim(res$reml$invAI), c(4L, 4L))

  # the reported estimate is the plug-in value implied by the variance
  # components, not the mean of the Monte Carlo draws (issue #13)
  v <- res$var[, 'Estimated variances']
  expect_equal(unname(res$funvars['mean', 'Heritability']),
               unname(v['genetic'] / sum(v)),
               tol = 1e-04)

  # heritability shown in summary, with all three reported numbers
  expect_output(print(summary(res)), "Heritability")
  expect_output(print(summary(res)), "Estimate")
  expect_output(print(summary(res)), "Sample Mean")

  # reported SE are consistent with AI matrix
  expect_equal(res$var[, 'S.E.'], sqrt(diag(res$reml$invAI)),
               tol = 1e-04, check.attributes = FALSE)
  
})


test_that('heritability and additional function are parsed correctly', {
  
  ## Simulate a dataset with a heritability of 2/(1+2+1+1) = 0.4
  set.seed(1234)
  dat <- breedR.sample.phenotype(fixed   = c(mu = 10, x = 2),
                                 genetic = list(model    = 'add_animal',
                                                Nparents = c(10, 10),
                                                sigma2_a = 2,
                                                check.factorial = FALSE),
                                 N = 1e3,
                                 residual.variance = 1)
  
  res <- remlf90(
    phenotype ~ 1 + X.x,
    genetic = list(model = 'add_animal',
                   pedigree = dat[, 1:3],
                   id = 'self'),
    progsf90.options = 'se_covar_function Halt G_3_3_1_1/(1+G_3_3_1_1+R_1_1)',
    data   = dat
  )
  
  expect_true(any(grepl('* SE for function of \\(co\\)variances Heritability', 
                        res$reml$output)))
  expect_true(any(grepl('Heritability  - Function: ', res$reml$output)))
  expect_true(any(grepl('Halt  - Function: ', res$reml$output)))

  expect_is(res$funvars, 'matrix')
  expect_is(res$reml$invAI, 'matrix')
  
  expect_identical(dim(res$reml$invAI), c(2L, 2L))
  
  expect_output(print(summary(res)), 'Halt')
  expect_output(print(summary(res)), 'Heritability')
})


test_that('AI-remlf90() without genetic does not return heritability but does return inverse AI matrix', {
  
  ## Simulate a small dataset for testing purposes
  dat <- breedR.sample.phenotype(fixed   = c(mu = 10),
                                 N = 100,
                                 residual.variance = 1)
  
  res <- remlf90(phenotype ~ 1,
                 data   = dat)
  
  expect_false(any(grepl('* SE for function of \\(co\\)variances', 
                         res$reml$output)))
  expect_false(any(grepl('  - Function: ', res$reml$output)))
  
  expect_identical(res$funvars, list())
  
  expect_is(res$reml$invAI, 'matrix')
  expect_identical(dim(res$reml$invAI), c(1L, 1L))
  
  expect_output(print(summary(res)), 'Variance components')
})


test_that('EM-remlf90() returns empty heritability and no inverse AI matrix', {

  ## Simulate a small dataset for testing purposes
  dat <- breedR.sample.phenotype(fixed   = c(mu = 10),
                                 genetic = list(model    = 'add_animal',
                                                Nparents = c(10, 10),
                                                sigma2_a = 2,
                                                check.factorial = FALSE),
                                 N = 100,
                                 residual.variance = 1)
  
  res <- remlf90(phenotype ~ 1 ,
                 genetic = list(model = 'add_animal',
                                pedigree = dat[, 1:3],
                                id = 'self'),
                 data   = dat,
                 method = 'em')
  
  expect_false(any(grepl('* SE for function of \\(co\\)variances Heritability', 
                         res$reml$output)))
  expect_false(any(grepl('Heritability  - Function: ', res$reml$output)))

  expect_identical(res$funvars, list())
  expect_null(res$reml$invAI)
  
  expect_output(print(summary(res)), 'Variance components:')
})


test_that('a genetic model fits with heterogeneous residual variances (issue #1)', {

  ## Data file columns: phe_X, Intercept, xc, genetic, so the covariate of
  ## the residual variance is column 3. x is rescaled to [0, 0.93] so that a
  ## slope of order 1 is a sensible start. A zero initial slope crashes
  ## BLUPF90+ 2.73, hence 0.1.
  dat <- globulus
  dat$xc <- dat$x / 100
  gen <- list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self')

  expect_error(
    res <- suppressMessages(
      remlf90(phe_X ~ xc, genetic = gen, data = dat,
              progsf90.options = hetres_options(covariate_cols = 3,
                                                initial = c(log(5), 0.1)))
    ),
    NA
  )

  ## the default heritability, which BLUPF90+ refuses under hetres, is not
  ## requested
  parameters_file <- readLines(file.path(res$reml$dir, 'parameters'))
  expect_false(any(grepl('se_covar_function', parameters_file)))
  expect_true(any(grepl('^OPTION hetres_pos 3', parameters_file)))
  expect_length(res$funvars, 0L)

  ## the coefficients replace the Residual row
  expect_identical(rownames(res$var), 'genetic')
  expect_identical(rownames(res$hetres), c('a0', 'a1'))
  expect_true(all(is.finite(res$hetres[, 'Estimate'])))
  expect_true(all(res$hetres[, 'S.E.'] > 0))
  expect_identical(rownames(res$reml$invAI), c('genetic', 'a0', 'a1'))
  expect_equal(unname(res$hetres[, 'S.E.']),
               unname(sqrt(diag(res$reml$invAI))[2:3]))

  expect_output(print(summary(res)), 'Residual variance model')

  ## control: without hetres the same model keeps its heritability and has
  ## no hetres element
  ctl <- suppressMessages(remlf90(phe_X ~ xc, genetic = gen, data = dat))
  expect_identical(colnames(ctl$funvars), 'Heritability')
  expect_null(ctl$hetres)
  expect_identical(rownames(ctl$var), c('genetic', 'Residual'))
})


test_that('heterogeneous residual variance coefficients are recovered', {

  ## log(var(e)) = 0.5 + 0.8 x, with a random group effect. At this seed the
  ## estimates sit at about 1.2 and 0.4 S.E. from the truth.
  set.seed(1)
  n <- 2000
  ng <- 50
  dat <- data.frame(g = factor(sample(ng, n, replace = TRUE)),
                    x = runif(n, 0, 2))
  u <- rnorm(ng, sd = sqrt(2))
  dat$y <- 10 + dat$x + u[dat$g] + rnorm(n, sd = sqrt(exp(0.5 + 0.8 * dat$x)))

  ## data file columns: y, Intercept, x, g
  res <- suppressMessages(
    remlf90(y ~ x, random = ~ g, data = dat,
            progsf90.options = hetres_options(covariate_cols = 3,
                                              initial = c(log(4), 0.1)))
  )

  z <- (res$hetres[, 'Estimate'] - c(a0 = 0.5, a1 = 0.8)) / res$hetres[, 'S.E.']
  expect_true(all(abs(z) < 3), label = paste('z =', toString(round(z, 2))))
  expect_identical(rownames(res$var), 'g')
})
