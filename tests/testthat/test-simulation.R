
#### Context: Simulation infrastructure ####
context("Simulation infrastructure") 

dat <- try(
  breedR.sample.phenotype(
    fixed   = c(mu = 10, x = 2),
    random = list(u = list(nlevels = 3,
                           sigma2  = 1)),
    genetic = list(model    = 'competition',
                   Nparents = c(10, 10),
                   sigma2_a = matrix(c(2, -1, -1, 2), 2, 2),
                   competition_decay = 1,
                   check.factorial = FALSE,
                   pec = 0.5),
    spatial = list(model     = 'AR',
                   grid.size = c(10, 10),
                   rho   = c(.9, .5),
                   sigma2_s  = 1),
    residual.variance = 1)
)

test_that('breedR.sample.phenotype() runs without error', {
  expect_false(inherits(dat, 'try-error'))
})

test_that("breedR.sample.phenotype() drops unsampled founders without a spatial effect (#61)", {
  ## 3 offspring from 4 dams and 4 sires: at most 3 of each are sampled, so
  ## whatever the seed, some founders have no offspring and must be removed
  N  <- 3
  Np <- c(4, 4)
  set.seed(61)
  dat <- suppressWarnings(breedR.sample.phenotype(
    fixed   = c(mu = 10, x = 2),
    random  = list(u = list(nlevels = 2, sigma2 = 1)),
    genetic = list(model    = 'add_animal',
                   Nparents = Np,
                   sigma2_a = 1,
                   check.factorial = FALSE),
    N = N))

  ## Replay the draws made before the pedigree, and the pedigree itself, to
  ## recover which simulated individual each returned row should be
  Nfull <- N + sum(Np)
  set.seed(61)
  x    <- stats::runif(Nfull)
  lev  <- sample(2, Nfull, replace = TRUE)
  val  <- stats::rnorm(2)
  dad  <- sample(Np[2], N, replace = TRUE)
  mum  <- sample(Np[1], N, replace = TRUE) + Np[2]
  kept <- sort(unique(c(dad, mum)))          # founders with offspring
  orig <- c(kept, sum(Np) + seq_len(N))      # original code of each row
  off  <- length(kept) + seq_len(N)          # offspring rows
  u    <- as.numeric(levels(dat$u))[dat$u]

  expect_lt(length(kept), sum(Np))
  expect_identical(names(dat), c('self', 'sire', 'dam', 'X.mu', 'X.x', 'u',
                                 'BV', 'resid', 'phenotype'))
  expect_equal(nrow(dat), length(orig))
  expect_equal(dat$self, seq_len(nrow(dat)))
  expect_true(all(is.na(dat$sire[-off]) & is.na(dat$dam[-off])))
  expect_equal(orig[dat$sire[off]], dad)
  expect_equal(orig[dat$dam[off]], mum)
  expect_equal(dat$X.x, x[orig])
  expect_equal(u, val[lev][orig])
  expect_equal(dat$phenotype, 10 + 2 * dat$X.x + u + dat$BV + dat$resid)
})
