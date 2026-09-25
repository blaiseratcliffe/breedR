
data(globulus)
ped <- build_pedigree(1:3, data = globulus)
# Test function
fit.model <- function(vi, vigen, random, dat = globulus, ...) {
  try(
    suppressMessages(
      remlf90(fixed   = phe_X ~ gen,
              random  = random,
              var.ini = vi,
              genetic = list(model = 'add_animal', 
                             var.ini = vigen,
                             pedigree = ped,
                             id = 'self'), 
              data    = dat)
    ),
  silent = TRUE)
}

# Test data
testdat <- list(
  list(
    vi     = NULL,   # Missing specification: FAIL
    vigen  = 1,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1),   # Missing specification of residual: FAIL
    vigen  = 1,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = list(gg = 1,      # Missing specification of bloc: FAIL
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1,     # Missing specification of gg: FAIL
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1,     # Missing genetic specification: FAIL
                  resid = 1),
    vigen  = NULL,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = NULL,              # OK: no specification at all
    vigen  = NULL,
    random = ~ bl,
    expectation = 1
  ),
  list(
    vi     = list(bl = 1,       # OK
                  resid = 1),
    vigen  = 1,
    random = ~ bl,
    expectation = 1
  ),
  list(
    vi     = list(bl = 1,       # OK
                  gg = 1,
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 1
  ),
  list(
    vi     = list(resid = 1),     # OK
    vigen  = 1,
    random = NULL,
    expectation = 1
  )
)


#### Context: Variance components specifications ####
context("Variance components specifications")

# reml results
# fit.model(vi=list(resid = 1), vigen=1, random = NULL)
# do.call('fit.model', testdat[[1]])
# do.call('fit.model', testdat[[7]])
reslst <- lapply(testdat, function(x) do.call(fit.model, x))


# Compare expected and true results
run_expectations <- function(m, res) {
  # Check that remlf90 behaves as expected
  test_that("remlf90 requires either full or null variance specifications", {
    ifelse( m$expectation,
            expect_true(!inherits(res, "try-error")),
            expect_true(inherits(res, "try-error")) )
  })
}

for(i in seq_along(testdat)) {
#   cat(i)
  run_expectations(testdat[[i]], reslst[[i]])
}

#### Context: Multitrait specifications ####
context("Multitrait interface")

## two correlated variables
dim <- 3
Nobs <- 1e4
Nbl <- 50
beta_X <- c(-1, 5, 3)
sample_covar <- function(dim) {
  x <- sample(-2:dim, size = dim**2, replace = TRUE)
  crossprod(matrix(x, nrow = dim))
}
set.seed(123)
S_bl <- sample_covar(dim)   # 5 & 3 & 5 // 22 & 10 // 11
S_resid <- sample_covar(dim)  # 9 & 3 & -3 // 9 & 9 // 14
# diag(1/sqrt(diag(S_bl))) %*% S_bl %*% diag(1/sqrt(diag(S_bl)))

dimnames(S_bl) <- dimnames(S_resid) <- 
  rep(list(paste0("y", seq_len(dim))), 2)

bl_levels <- paste0(
  "bl",
  sprintf(paste0("%0", floor(log10(Nbl)+1), "d"), seq_len(Nbl))
)

testdat <- data.frame(
  X = runif(Nobs),
  breedR.sample.ranef(
    dim, S_bl, Nbl, labels = bl_levels, N = Nobs, vname = 'bl'
  ),
  breedR.sample.ranef(dim, S_resid, Nobs, vname = 'e'))

# var(testdat[, c('e1', 'e2')])  # ~ S_resid

testdat <- 
  transform(testdat,
            y1 = beta_X[1]*X + bl_y1 + e_y1,
            y2 = beta_X[2]*X + bl_y2 + e_y2,
            y3 = beta_X[3]*X + bl_y3 + e_y3)


test_that("Residual variance acurately identified in a fixed-effects model", {
  
  ## All fixed effects (AI fails with 3 traits)
  res <- remlf90(
    cbind(y1, y2, y3) ~ X + bl,
    data = testdat,
    method = "em"
  )
  
  expect_equal(S_resid, res$var$Residual, tol = .01)
})


test_that("Simulated values reasonably recovered using one or more traits", {
  
  ## 1 trait
  res_1 <- remlf90(
    cbind(y1) ~ 0 + X,
    random = ~ bl,
    data = testdat
  )
  
  expect_equal(S_resid[1, 1], res_1$var["Residual", 1], tol = .01)
  expect_equal(S_bl[1, 1], res_1$var["bl", 1], tol = .1)
  expect_equal(beta_X[1], fixef(res_1)$X, tol = .1, check.attributes = FALSE)

  ## 2 trait
  res_2 <- remlf90(
    cbind(y1, y2) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "ai"
  )

  expect_equal(S_resid[-3, -3], res_2$var[["Residual", 1]], tol = .01)
  expect_equal(S_bl[-3, -3], res_2$var[["bl", 1]], tol = 1)
  expect_equal(beta_X[-3], fixef(res_2)$X, tol = .1, check.attributes = FALSE)
  
  ## 3 trait  (with "em" since otherwise does not converge)
  res_3 <- remlf90(
    cbind(y1, y2, y3) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "em"
  )  

  expect_equal(S_resid, res_3$var$Residual, tol = .01)
  expect_equal(S_bl, res_3$var$bl, tol = 1)
  expect_equal(beta_X, fixef(res_3)$X, tol = .01, check.attributes = FALSE)
})


test_that("Initial variance specification", {
  
  vi <- function(n) diag(n) + matrix(1, n, n)
  ## 2 trait - full matrices
  res_2 <- remlf90(
    cbind(y1, y2) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "ai",
    var.ini = list(bl = vi(2), resid = vi(2))
  )
  
  expect_equal(S_resid[-3, -3], res_2$var[["Residual", 1]], tol = .01)
  expect_equal(S_bl[-3, -3], res_2$var[["bl", 1]], tol = 1)
  expect_equal(beta_X[-3], fixef(res_2)$X, tol = .1, check.attributes = FALSE)
  
  ## The invAI matrix have all covariance terms
  invai.names <- c(paste0("bl.", c("y1", "y1_bl.y2", "y2")),
                   paste0("resid.", c("y1", "y1_resid.y2", "y2")))
  expect_identical(rownames(res_2$reml$invAI), invai.names)

  ## 2 trait - non-full residual matrix
  res_2 <- remlf90(
    cbind(y1, y2) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "ai",
    var.ini = list(bl = vi(2), resid = diag(2))
  )
  
  expect_equal(diag(diag(S_resid[-3, -3])), res_2$var[["Residual", 1]],
               tol = .01, check.attributes = FALSE)
  
  ## estimated residual covariance of 0
  expect_identical(res_2$var[["Residual", 1]][1, 2], 0)
  
  expect_equal(S_bl[-3, -3], res_2$var[["bl", 1]], tol = 1)
  expect_equal(beta_X[-3], fixef(res_2)$X, tol = .1, check.attributes = FALSE)
  
  ## The invAI matrix does not have the residual covariance term
  invai.names <- c(paste0("bl.", c("y1", "y1_bl.y2", "y2")),
                   paste0("resid.", paste0("y", 1:2)))
  expect_identical(rownames(res_2$reml$invAI), invai.names)
  
  ## 2 trait - non-full block matrix
  res_2 <- remlf90(
    cbind(y1, y2) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "ai",
    var.ini = list(bl = diag(2), resid = vi(2))
  )
  
  expect_equal(S_resid[-3, -3], res_2$var[["Residual", 1]],
               tol = .01, check.attributes = FALSE)

  expect_equal(S_bl[-3, -3], res_2$var[["bl", 1]], tol = 1)
  
  ## estimated bl covariance of 0
  expect_identical(res_2$var[["bl", 1]][1, 2], 0)
  
  expect_equal(beta_X[-3], fixef(res_2)$X, tol = .1, check.attributes = FALSE)

  ## The invAI matrix does not have the block covariance term
  invai.names <- c(paste0("bl.", c("y1", "y2")),
                   paste0("resid.", c("y1", "y1_resid.y2", "y2")))
  expect_identical(rownames(res_2$reml$invAI), invai.names)

  
  ## 3 trait - not full matrices
  blvi <- vi(3)
  blvi[2, 1] <- blvi[1, 2] <- 0
  
  res_3 <- remlf90(
    cbind(y1, y2, y3) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "em",
    var.ini = list(bl = blvi, resid = diag(3))
  )  
  
  expect_equal(diag(diag(S_resid)), res_3$var$Residual,
               tol = .1, check.attributes = FALSE)
  expect_equal(S_bl, res_3$var$bl, tol = 1)

  ## estimated bl covariances of 0
  expect_identical(res_3$var$bl[1, 2], 0)
  
  ## estimated residual covariance of 0
  expect_identical(res_3$var$Residual[lower.tri(res_3$var$Residual)], rep(0, 3))

  expect_equal(beta_X, fixef(res_3)$X, tol = .5, check.attributes = FALSE)
})


# mf <- model.frame(cbind(V1, V2) ~ 0 + mu, transform(testdat, mu = 1))
# attr(attr(mf, 'terms'), 'term.types') <- list(mu = "fixed")
# eff <- build.effects(mf, genetic = NULL, spatial = NULL, generic = NULL, var.ini = S)
# pf90 <- progsf90(mf, eff, res.var.ini = S)

## Use larix dataset:
## - Two phenotypes: LAS and DOS
## - repeated measurements along 16 years (yr)


test_that("Multitrait model with all kind of effects works as expected", {
  
  inc.mat <- model.matrix(~ 0 + bl, larix)
  cov.mat <- diag(nlevels(larix$bl))
  
  fullrun <- function(method, opt = NULL) {
    try(
        remlf90(
          fixed   = cbind(LAS, DOS) ~ rep,
          random  = ~ bl,
          genetic = list(model = 'add_animal',
                         pedigree = larix[, 1:3],
                         id = 'self'),
          spatial = list(model = 'AR',
                         coordinates = larix[, c('x', 'y')],
                         rho = c(.8, .8)),
          generic = list(block = list(inc.mat,
                                      cov.mat)),
          data    = larix,
          method = method,
          progsf90.options = opt
        )
    )
  }
  
  ## make things fast, as I am not looking at numerical results
  res_ai <- fullrun("ai", opt = c("maxrounds 2"))
  ## cannot use it with em as the logfile would not report final estimates
  res_em <- fullrun("em")
  
  fixef_names <- "rep"
  ranef_names <- c("bl", "genetic", "spatial", "block")
  
  
  ## No errors
  expect_false(inherits(res_em, "try-error"))
  expect_false(inherits(res_ai, "try-error"))
  
  ## fixed effect estimates
  expect_identical(names(fixef(res_em)), fixef_names)
  expect_identical(names(fixef(res_ai)), fixef_names)
  
  ## variance component estimates
  ## em: a list -> names
  ## ai: a matrix (effects x (estimate, se)) -> rownames
  expect_identical(names(res_em$var), c(ranef_names, "Residual"))
  expect_identical(rownames(res_ai$var), c(ranef_names, "Residual"))
  
  ## random effect blups
  expect_identical(names(ranef(res_em)), ranef_names)
  expect_identical(names(ranef(res_ai)), ranef_names)
  
})


test_that("trait-specific random effects match balanced REML and BLUPs", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")
  a <- c(-3, -2, -1, 1, 2, 3)
  dat <- data.frame(
    rep = factor(rep(seq_along(a), each = 4)),
    y1 = 10 + rep(a, each = 4) + rep(c(-1, 1, -1, 1), 6),
    y2 = 20 + rep(c(-2, -2, 2, 2), 6)
  )
  opts <- c("maxrounds 100", "conv_crit 1d-12")

  ## The residual covariance is constrained to zero, so the joint likelihood
  ## factors into a balanced random-intercept fit and an intercept-only fit.
  ## MSwithin = 4/3, MSbetween = 112/5, G = (MSbetween - MSwithin)/4.
  for (method in c("ai", "em")) {
    fit <- suppressMessages(remlf90(
      cbind(y1, y2) ~ 1, random = ~ rep, data = dat,
      traits = list(rep = "y1"), method = method,
      var.ini = list(rep = diag(c(5, 0)), residuals = diag(c(1, 4))),
      progsf90.options = opts
    ))
    g <- if (method == "ai") fit$var[["rep", 1]] else fit$var$rep
    r <- if (method == "ai") fit$var[["Residual", 1]] else fit$var$Residual
    expect_equal(unname(g[1, 1]), 79/15, tolerance = 1e-4)
    expect_true(all(is.na(g[2, ])))
    expect_true(all(is.na(g[, 2])))
    expect_equal(unname(diag(r)), c(4/3, 96/23), tolerance = 1e-4)
    expect_identical(unname(r[1, 2]), 0)
    expect_equal(as.numeric(fixef(fit)$Intercept), c(10, 20),
                 tolerance = 1e-5)
    expect_equal(as.numeric(ranef(fit)$rep[, "y1"]), 79/84*a,
                 tolerance = 1e-4)
    expect_true(all(is.na(ranef(fit)$rep[, "y2"])))
    expect_output(print(summary(fit)), 'Absent random effects: rep on y2')
    expect_equal(as.numeric(fitted(fit)[, 1]),
                 10 + rep(79/84*a, each = 4), tolerance = 1e-4)
    expect_equal(as.numeric(fitted(fit)[, 2]), rep(20, 24),
                 tolerance = 1e-5)

    first <- suppressMessages(remlf90(
      y1 ~ 1, random = ~ rep, data = dat, method = method,
      var.ini = list(rep = 5, residuals = 1), progsf90.options = opts
    ))
    second <- suppressMessages(remlf90(
      y2 ~ 1, data = dat, method = method,
      var.ini = list(residuals = 4), progsf90.options = opts
    ))
    expect_equal(as.numeric(fitted(fit)[, 1]), as.numeric(fitted(first)),
                 tolerance = 1e-4)
    expect_equal(as.numeric(fitted(fit)[, 2]), as.numeric(fitted(second)),
                 tolerance = 1e-5)
    if (method == "ai") {
      expect_identical(rownames(fit$reml$invAI),
                       c("rep.y1", "resid.y1", "resid.y2"))
      expect_true(all(is.na(fit$var[["rep", 2]][2, ])))
    }
  }
})


test_that("default and all-trait selections preserve parameters and estimates", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")
  a <- rep(c(-3, -2, -1, 1, 2, 3), each = 4)
  dat <- data.frame(rep = factor(rep(1:6, each = 4)),
                    y1 = 10 + a + rep(c(-1, 1, -1, 1), 6),
                    y2 = 20 + .75*a + rep(c(-2, -2, 2, 2), 6))
  fit_default <- function(...) suppressMessages(remlf90(
    cbind(y1, y2) ~ 1, random = ~ rep, data = dat,
    var.ini = list(rep = diag(c(5, 2)), residuals = diag(c(1, 4))),
    progsf90.options = c("maxrounds 100", "conv_crit 1d-12"), ...
  ))
  original <- fit_default()
  parameters <- readLines(file.path(original$reml$dir, "parameters"))
  for (selection in list(NULL, list(), list(rep = c("y2", "y1")))) {
    fit <- fit_default(traits = selection)
    expect_identical(readLines(file.path(fit$reml$dir, "parameters")), parameters)
    expect_identical(fit$effects, original$effects)
    expect_equal(fit$var, original$var, tolerance = 1e-12)
    expect_equal(fitted(fit), fitted(original), tolerance = 1e-12)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(original)),
                 tolerance = 1e-12)
  }
})


test_that("noncontiguous traits agree with a hand-written BLUPF90 model", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")
  set.seed(51)
  group <- rep(seq_len(20), each = 6)
  u1 <- rnorm(20, sd = 1.5)
  u3 <- .4*u1 + rnorm(20, sd = 1.3)
  shared <- rnorm(length(group), sd = .5)
  dat <- data.frame(
    y1 = 10 + u1[group] + shared + rnorm(length(group), sd = 1.5),
    y2 = 20 + shared + rnorm(length(group), sd = 1.5),
    y3 = 30 + u3[group] + shared + rnorm(length(group), sd = 1.5),
    rep = factor(group)
  )
  initial_g <- matrix(c(2, 0, .4, 0, 0, 0, .4, 0, 2), 3)
  initial_r <- matrix(.5, 3, 3) + diag(2.5, 3)
  fit <- suppressMessages(remlf90(
    cbind(y1, y2, y3) ~ 1, random = ~ rep, data = dat,
    traits = list(rep = c("y3", "y1")),
    var.ini = list(rep = initial_g, residuals = initial_r),
    progsf90.options = c("maxrounds 100", "conv_crit 1d-10")
  ))

  ## Independent serialization: no breedR render or parse helper is used here.
  direct_dir <- tempfile("trait-presence-direct-")
  dir.create(direct_dir)
  on.exit(unlink(direct_dir, recursive = TRUE), add = TRUE)
  write.table(data.frame(dat[, 1:3], intercept = 1, group = group),
              file.path(direct_dir, "data"), row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  writeLines(c(
    "DATAFILE", "data", "NUMBER_OF_TRAITS", "3", "NUMBER_OF_EFFECTS", "2",
    "OBSERVATION(S)", "1 2 3", "WEIGHT(S)", "", "EFFECTS:",
    "4 4 4 1 cross", "5 0 5 20 cross", "RANDOM_RESIDUAL VALUES",
    "3 .5 .5", ".5 3 .5", ".5 .5 3", "RANDOM_GROUP", "2",
    "RANDOM_TYPE", "diagonal", "FILE", "", "(CO)VARIANCES",
    "2 0 .4", "0 0 0", ".4 0 2", "OPTION method VCE",
    "OPTION maxrounds 100", "OPTION conv_crit 1d-10", "OPTION sol se"
  ), file.path(direct_dir, "parameters"))
  binary <- file.path(breedR.getOption("breedR.bin"),
                      breedR:::progsf90_files(breedR:::breedR.os.type()))
  oldwd <- setwd(direct_dir)
  direct <- tryCatch(system2(binary, input = "parameters", stdout = TRUE),
                     finally = setwd(oldwd))
  expect_true(file.exists(file.path(direct_dir, "solutions")))
  raw_solutions <- read.table(file.path(direct_dir, "solutions"), skip = 1)
  for (trait in c(1, 3)) {
    rows <- raw_solutions[raw_solutions$V1 == trait & raw_solutions$V2 == 2, ]
    rows <- rows[order(rows$V3), ]
    expect_equal(as.numeric(ranef(fit)$rep[, trait]), rows$V4, tolerance = 1e-4)
  }
  read_block <- function(label) {
    start <- tail(grep(label, direct, fixed = TRUE), 1)
    unname(as.matrix(read.table(text = paste(direct[start + 1:3], collapse = "\n"))))
  }
  expect_equal(unname(fit$var[["rep", 1]][c(1, 3), c(1, 3)]),
               read_block("Genetic variance(s)")[c(1, 3), c(1, 3)],
               tolerance = 1e-4)
  expect_equal(unname(fit$var[["Residual", 1]]), read_block("Residual variance(s)"),
               tolerance = 1e-4)
  expect_true(all(is.na(ranef(fit)$rep[, "y2"])))
  expect_identical(dim(fit$reml$invAI), c(9L, 9L))
})


test_that("site-as-trait fit without var.ini (#71)", {
  skip_if_not(isTRUE(check_progsf90(quiet = TRUE)), "PROGSF90 binaries not installed")

  ## 30 sires with 5 offspring at each of two sites; every record is observed
  ## at one site only, so no record observes both traits
  set.seed(71)
  ns <- 30; noff <- 5
  a_sire <- matrix(rnorm(2 * ns), ns) %*% chol(matrix(c(1, .7, .7, 1), 2))
  ped <- data.frame(self = 1:(ns + 2 * ns * noff),
                    sire = c(rep(0, ns), rep(rep(1:ns, each = noff), 2)),
                    dam = 0)
  dat <- data.frame(self = (ns + 1):(ns + 2 * ns * noff),
                    sire = rep(rep(1:ns, each = noff), 2),
                    site = rep(1:2, each = ns * noff))
  bv <- 0.5 * a_sire[cbind(dat$sire, dat$site)] +
    rnorm(nrow(dat), sd = sqrt(0.75))
  dat$y1 <- ifelse(dat$site == 1, 10 + bv + rnorm(nrow(dat), sd = sqrt(2.3)), NA)
  dat$y2 <- ifelse(dat$site == 2, 20 + bv + rnorm(nrow(dat), sd = sqrt(2.3)), NA)
  expect_false(any(complete.cases(dat[, c("y1", "y2")])))

  expect_message(
    fit <- remlf90(cbind(y1, y2) ~ 1,
                   genetic = list(model = "add_animal", pedigree = ped,
                                  id = "self"),
                   data = dat, method = "ai"),
    "Using default initial variances")
  G <- fit$var[["genetic", "Estimated variances"]]
  R <- fit$var[["Residual", "Estimated variances"]]
  expect_true(all(is.finite(G)))
  expect_true(all(eigen(G, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_true(G[1, 2] != 0)
  ## the residual covariance of the two sites cannot be estimated
  expect_identical(unname(R[1, 2]), 0)

  ## The same starting values given explicitly, with the residual covariance
  ## at 0: the same fit, and no extra parameter in logLik()'s df
  G0 <- breedR:::default_initial_variance(dat[, c("y1", "y2")],
                                          cor.effect = 0.1, digits = 2)
  explicit <- suppressMessages(
    remlf90(cbind(y1, y2) ~ 1,
            genetic = list(model = "add_animal", pedigree = ped, id = "self",
                           var.ini = G0),
            var.ini = list(residuals = diag(diag(G0))),
            data = dat, method = "ai"))
  expect_identical(attr(logLik(fit), "df"), attr(logLik(explicit), "df"))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(explicit)))
})
