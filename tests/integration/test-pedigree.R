#### pedigree building and checking ####

context("Pedigree")

# Toy dataset with silly pedigree
test.dat <- data.frame(matrix(sample(100, 15), 5, 3,
                              dimnames = list(NULL, c('self', 'sire', 'dam'))),
                       y = rnorm(5))
ped.fix <- suppressWarnings(build_pedigree(1:3, data = test.dat))
test.res <- try(
  suppressMessages(
    suppressWarnings(
      remlf90(y~1,
              genetic = list(model = 'add_animal',
                             pedigree = test.dat[, 1:3],
                             id = 'self'),
              data = test.dat)
    )
  ),
  silent = TRUE
)

test_that('remlf90() builds and recodes the pedigree', {
  expect_false(inherits(test.res, 'try-error'))
})

test_that('get_pedigree() returns the recoded pedigree', {
  expect_identical(ped.fix, get_pedigree(test.res))
})


# Check that remlf90 handles correctly recoded pedigrees
# by comparing the genetics evaluations of a dataset with or without
# a shuffled pedigree

data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

res_ok <- try(
  suppressMessages(
    remlf90(fixed = phe_X ~ sex, 
            genetic = list(model = 'add_animal', 
                           pedigree = ped,
                           id = 'self'), 
            data = dat)
  )
)

# Shuffle the pedigree
mcode <- max(as.data.frame(ped), na.rm = TRUE)
map <- rep(NA, mcode)
set.seed(1234)
map <- sample(10*mcode, size = mcode)
m1_shuffled <- m1
m1_shuffled$Data[, 1:3] <- sapply(as.data.frame(ped), function(x) map[x])

ped_fix <- suppressWarnings(
  build_pedigree(1:3, data = as.data.frame(get_pedigree(m1_shuffled)))
)


res_shuffled <- try(
  suppressMessages(
    remlf90(fixed = phe_X ~ sex,
            genetic = list(model = 'add_animal', 
                           pedigree = ped_fix,
                           id = 'self'), 
            data = as.data.frame(m1_shuffled))
  )
)

# Except the call, and the reml output everything must be the same
# Update: also need to omit the shuffled random effects estimations
# which should be the same, but reordered
test_that('remlf90 handles recoded pedigrees correctly', {
  omit.idx <- match(c('call', 'effects', 'reml', 'ranef'), names(res_ok))
  expect_that(res_ok[-omit.idx], equals(res_shuffled[-omit.idx]))
})

test_that("gibbsf90() returns the pedigree that translates its genetic levels (#64)", {
  ## reversed codes put offspring before parents: the pedigree is recoded
  N <- max(globulus[, c('self', 'dad', 'mum')])
  flip <- function(x) ifelse(x == 0, 0L, as.integer(N + 1L - x))
  gr <- globulus
  gr[, c('self', 'dad', 'mum')] <-
    lapply(globulus[, c('self', 'dad', 'mum')], flip)
  dat <- gr[1:300, ]
  gen <- list(model = 'add_animal', pedigree = gr[, 1:3], id = 'self')
  gb <- suppressWarnings(suppressMessages(
    gibbsf90(phe_X ~ 1, genetic = gen, data = dat, n_samples = 2000L,
             burnin = 500L, thin = 10L, seed = c(11, 22))))
  on.exit(unlink(gb$dir, recursive = TRUE), add = TRUE)
  ped <- suppressWarnings(build_pedigree(1:3, data = gr[, 1:3]))
  expect_false(is.null(attr(ped, 'map')))
  expect_identical(gb$pedigree, ped)

  ## level k of the genetic effect is the animal coded k
  sol <- gb$solutions[gb$solutions$effect == 2, ]
  expect_identical(as.integer(sol$level), seq_along(ped@label))
  id <- match(sol$level, attr(gb$pedigree, 'map'))

  ## read through the map, the posterior means belong to the animals of
  ## the REML BLUPs: cor 0.98 here, 0.01 with the level taken as the id
  reml <- suppressWarnings(suppressMessages(
    remlf90(phe_X ~ 1, genetic = gen, data = dat)))
  bv <- ranef(reml)$genetic
  expect_gt(cor(sol$solution[match(as.integer(names(bv)), id)],
                as.numeric(bv)), 0.9)
})
