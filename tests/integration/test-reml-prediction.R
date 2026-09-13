### For testing prediction, we perform a cross-validation excercise ###

data(m1)

# Sample individuals and generate a secondary response
# with removed observations
Nobs <- nrow(as.data.frame(m1))
sel.idx <- sample(as.data.frame(m1)$self, Nobs/6)

m1$Data$y <- m1$Data$phe_X
m1$Data$y[sel.idx] <- NA
dat <- as.data.frame(m1)

# Fit with no missings
res.full <- try(
  suppressMessages(
    remlf90(fixed   = phe_X ~ sex, 
            genetic = list(model = 'add_animal', 
                           pedigree = get_pedigree(m1),
                           id = 'self'), 
            spatial = list(model = 'AR', 
                           coord = coordinates(m1),
                           rho = c(.9, .9)), 
            data = dat,
            method = 'ai')
  )
)

# Fit with missings
res.pred <- try(
  suppressMessages(
    remlf90(fixed   = y ~ sex, 
            genetic = list(model = 'add_animal', 
                           pedigree = get_pedigree(m1),
                           id = 'self'), 
            spatial = list(model = 'AR', 
                           coord = coordinates(m1),
                           rho = c(.9, .9)), 
            data = dat,
            method = 'ai')
  )
)

#### Context: Prediction and cross-validation ####
context("Prediction")

# For a cross-validation, the results will not be exactly equal
# but similar up to a given tolerance
tol = 2e-01

# Variance components
test_that("Estimated variance components are similar", {
  expect_equal(res.full$var, res.pred$var, tolerance = tol)
})

# Prediction of the spatial effect
# qplot(ranef(res.full)$spatial,
#       ranef(res.pred)$spatial) +
#   geom_abline(intercept=0, slope=1)
test_that("Predicted spatial effects are similar everywhere", {
  expect_equal(ranef(res.full)$spatial,
               ranef(res.pred)$spatial,
               tolerance = tol, check.attributes = FALSE)
})



# Prediction of the genetic effect
# plotdat <- data.frame(fullBV = ranef(res.full)$genetic,
#                       predBV = ranef(res.pred)$genetic,
#                       miss   = factor(0, levels = 0:1))
# plotdat$miss[sel.idx] <- 1L
# qplot(fullBV, predBV, data = plotdat, color = miss) + 
#   geom_abline(intercept = 0, slope = 1, col ='darkgray')
test_that("Predicted Breeding Values are similar in observed individuals", {
  expect_equal(ranef(res.full)$genetic[-sel.idx],
            ranef(res.pred)$genetic[-sel.idx],
            tolerance = tol)
})

test_that("Predicted Breeding Values of missings are still similar 
          up to one additional order of magnitude", {
  expect_equal(ranef(res.full)$genetic[sel.idx],
               ranef(res.pred)$genetic[sel.idx],
               tolerance = 10*tol)
})

# Estimation of fixed effects: can change up to one additional order of magnitude
test_that("Fixed effects are similar", {
  expect_equal(fixef(res.full), fixef(res.pred), tolerance = 10*tol)
})



test_that('(ai)remlf90() predict correctly when missing code is not 0', {
  ## dataset with positive and negative values
  dat <- breedR.sample.phenotype(fixed = c(mu = 0), N = 1e3)
  dat$group <- factor(rep(letters[1:4], each = 1e3/4))
  dat$phenotype <- dat$phenotype + as.numeric(dat$group)
  dat$phenotype[1] <- NA
  
  expect_error(
    res <- remlf90(phenotype ~ group, data = dat),
    NA
  )
  
  expect_equal(fitted(res)[1], fixef(res)$group[1], 
               check.attributes = FALSE)
})


test_that("trait absence preserves missing-response and factor-level alignment", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")
  a <- c(-3, -2, -1, 1, 2, 3)
  dat <- data.frame(
    rep = factor(rep(letters[1:6], each = 4), levels = rev(letters[1:6])),
    y1 = 10 + rep(a, each = 4) + rep(c(-1, 1, -1, 1), 6),
    y2 = 20 + rep(c(-2, -2, 2, 2), 6)
  )
  dat$y1[c(2, 11)] <- NA
  dat$y2[c(3, 17)] <- NA
  rownames(dat) <- paste0("tree", seq_len(nrow(dat)))
  fit_selected <- function(d) suppressMessages(remlf90(
    cbind(y1, y2) ~ 1, random = ~ rep, data = d,
    traits = list(rep = "y1"),
    var.ini = list(rep = diag(c(5, 0)), residuals = diag(c(1, 4))),
    progsf90.options = c("maxrounds 100", "conv_crit 1d-12")
  ))
  fit <- fit_selected(dat)
  expect_identical(rownames(ranef(fit)$rep), levels(dat$rep))
  expect_identical(rownames(model.frame(fit)), rownames(dat))
  expect_identical(dim(fitted(fit)), c(24L, 2L))
  expect_identical(rownames(fitted(fit)), rownames(dat))
  expect_identical(colnames(fitted(fit)), c("y1", "y2"))
  expect_true(all(is.finite(fitted(fit))))
  expect_equal(as.numeric(fitted(fit)[, 2]), rep(mean(dat$y2, na.rm = TRUE), 24),
               tolerance = 1e-5)
  expect_equal(unname(residuals(fit)),
               unname(as.matrix(dat[, c("y1", "y2")]) - fitted(fit)),
               check.attributes = FALSE)
  expect_identical(is.na(residuals(fit)),
                   is.na(as.matrix(dat[, c("y1", "y2")])))

  permutation <- c(seq(24, 2, by = -2), seq(23, 1, by = -2))
  shuffled <- fit_selected(dat[permutation, ])
  expect_equal(unname(fitted(shuffled)[order(permutation), ]), unname(fitted(fit)),
               tolerance = 1e-4)
  expect_equal(ranef(shuffled)$rep, ranef(fit)$rep, tolerance = 1e-4)
})


test_that("a restricted single-level generic effect preserves trait matrices", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")
  x <- seq(-2, 2, length.out = 24)
  dat <- data.frame(x = x,
                    y1 = 10 + 2*x + rep(c(-.1, .1), 12),
                    y2 = 3*x + rep(c(-1, 1), 12))
  fit <- suppressMessages(remlf90(
    cbind(y1, y2) ~ 0 + x, data = dat,
    generic = list(shared = list(incidence = matrix(1, 24, 1),
                                  covariance = matrix(1),
                                  var.ini = diag(c(100, 0)))),
    traits = list(shared = "y1"),
    var.ini = list(residuals = diag(c(.01, 1))),
    progsf90.options = c("maxrounds 100", "conv_crit 1d-10")
  ))
  shared <- ranef(fit)$shared
  expect_identical(dim(shared), c(1L, 2L))
  expect_identical(dim(attr(shared, "se")), c(1L, 2L))
  expect_identical(colnames(shared), c("y1", "y2"))
  expect_identical(dimnames(attr(shared, "se")), dimnames(shared))
  expect_true(is.finite(shared[1, "y1"]))
  expect_true(is.finite(attr(shared, "se")[1, "y1"]))
  expect_true(is.na(shared[1, "y2"]))
  expect_true(is.na(attr(shared, "se")[1, "y2"]))
  expect_identical(dim(fixef(fit)$x), c(1L, 2L))
  expect_identical(dim(fitted(fit)), c(24L, 2L))
  expect_true(all(is.finite(fitted(fit))))
  ## With no shared random effect on y2 and zero residual covariance, its
  ## prediction is the ordinary least-squares line through the origin.
  expected_y2 <- x * sum(x * dat$y2) / sum(x^2)
  expect_equal(as.numeric(fitted(fit)[, "y2"]), expected_y2, tolerance = 1e-5)
})
