#### pedigree building and checking ####

context("Pedigree infrastructure")

# Retrieve pedigree from remlf90 objects

test_that('get_pedigree() returns NULL when there is no genetic effect', {
  res <- load_res("fixonly")
  expect_true(is.null(get_pedigree(res)))
})


# Use the pedigree in data(m4) and shuffle the codes
data(m4)
ped <- as.data.frame(m4)[, c('self', 'dad', 'mum')]

test_that('The pedigree from m4 is not complete, but otherwise correct', {
  expect_true(!check_pedigree(ped)['full_ped'])
  expect_true(all(check_pedigree(ped)[-1]))
})

# Generate a crazy map
mcode <- max(ped, na.rm = TRUE)
map <- rep(NA, mcode)
set.seed(1234)
map <- sample(10*mcode, size = mcode)

# Generate a crazy pedigree that fails all checks
ped_shuffled <- sapply(ped, function(x) map[x])
# Introduce some unknown parents either with NA or with 0
ped_shuffled[, 2:3][sample(2*nrow(ped), 200)] <- c(0, NA)

test_that('The shuffled pedigree fails all checks', {
  expect_true(all(!check_pedigree(ped_shuffled)))
})

# Reorder and recode 
ped_fix <- suppressWarnings(build_pedigree(1:3, data = ped_shuffled))
test_that('build_pedigree() fixes everything', {
  expect_true(all(check_pedigree(ped_fix)))
})


test_that('build_pedigree() labels double codes of 1e5 and above in full (#43)', {

  ## Codes 1..n need no recoding, so the labels are the codes themselves.
  ## pedigreemm stores labels with as.character(), which writes a double
  ## 100000 as "1e+05", and then no integer id matches it.
  n <- 100001
  ped_dbl <- data.frame(self = as.numeric(seq_len(n)), dad = 0, mum = 0)
  ped_dbl[n, c('dad', 'mum')] <- c(1e5, 5e4)
  ped_int <- data.frame(self = seq_len(n), dad = 0L, mum = 0L)
  ped_int[n, c('dad', 'mum')] <- c(100000L, 50000L)

  ped43 <- build_pedigree(1:3, data = ped_dbl)

  expect_null(attr(ped43, 'map'))
  expect_identical(ped43@label[1e5], '100000')

  ## parents still point at the right individuals
  expect_identical(ped43@sire[n], 100000L)
  expect_identical(ped43@dam[n], 50000L)

  ## the storage type of the codes makes no difference
  expect_identical(ped43, build_pedigree(1:3, data = ped_int))
})

