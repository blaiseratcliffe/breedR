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

test_that('The pedigree from m4 is not complete, so its codes do not start at 1', {
  ## The founders are missing, so the lowest code is 161 (#50)
  checks <- check_pedigree(ped)
  expect_false(checks[['full_ped']])
  expect_false(checks[['codes_consecutive']])
  expect_true(all(checks[c('offsp_follows', 'codes_sorted')]))
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


test_that('check_pedigree() requires the codes to start at 1 (#50)', {
  expect_false(
    check_pedigree(data.frame(self = 2:4, dad = 0, mum = 0))[['codes_consecutive']])
  expect_true(
    check_pedigree(data.frame(self = 1:3, dad = 0, mum = 0))[['codes_consecutive']])
})


test_that('build_pedigree() recodes a pedigree coded from 2 (#50)', {

  ## Sorted, consecutive, and offspring follow parents. But the codes are used
  ## as positions 1..n, so a pedigree that starts at 2 must be recoded.
  ped2 <- data.frame(self = 2:5, dad = c(0L, 0L, 2L, 2L), mum = c(0L, 0L, 3L, 4L))
  expect_warning(ped50 <- build_pedigree(1:3, data = ped2), 'recoded')

  map <- attr(ped50, 'map')
  expect_identical(map, c(NA, 1:4))
  expect_true(all(check_pedigree(ped50)))

  ## translated back through the map, it is the pedigree we started from
  back <- as.data.frame(lapply(as.data.frame(ped50),
                               function(x) ifelse(is.na(x), 0L, match(x, map))))
  back <- back[order(back$self), ]
  rownames(back) <- NULL
  expect_identical(back, setNames(ped2, c('self', 'sire', 'dam')))
})


test_that('build_pedigree() refuses individual codes below 1 (#50)', {

  ## 0 marks an unknown parent, and the recoding map is indexed by code. An
  ## animal coded 0 used to go unrecoded, and every animal's data then went to
  ## the animal coded one below it.
  ped0 <- data.frame(self = 0:4, dad = c(0, 0, 0, 1, 1), mum = c(0, 0, 0, 2, 2))
  expect_error(build_pedigree(1:3, data = ped0), 'codes must be 1 or greater')

  ped_neg <- data.frame(self = -1:3, dad = c(0, 0, 0, 1, 1), mum = c(0, 0, 0, 2, 2))
  expect_error(build_pedigree(1:3, data = ped_neg), 'codes must be 1 or greater')
})


test_that('build_pedigree() refuses negative parent codes, and names them (#66)', {

  ## 0 or NA marks an unknown parent. A negative code is neither: with
  ## consecutive codes pedigreemm refused it with an unrelated message, and
  ## when the pedigree was recoded it became a negative subscript of the map.
  msg <- 'Negative parent codes: -1$'
  expect_error(build_pedigree(1:3, data = data.frame(self = 1:4,
                                                     dad = c(0, -1, 1, 1),
                                                     mum = c(0, 0, 2, 2))), msg)
  expect_error(build_pedigree(1:3, data = data.frame(self = c(1, 2, 3, 5),
                                                     dad = c(0, -1, 1, 1),
                                                     mum = c(0, 0, 2, 2))), msg)

  ## A negative subscript drops map entries. Dropping as many as the map has
  ## gaps gave a sire column of the right length, and animal 7 was silently
  ## given animal 5 as its sire.
  expect_error(build_pedigree(1:3, data = data.frame(self = c(1, 5, 6, 7),
                                                     dad = c(-1, -6, -7, -1),
                                                     mum = NA_real_)),
               'Negative parent codes: -7, -6, -1$')
})


test_that('build_pedigree() refuses codes above the integer range (#67)', {

  ## Codes are stored as integers, and the recoding map is indexed by code:
  ## a 10-digit id would need a map of billions of entries. Stop before that.
  big <- data.frame(self = 3e9 + 1:4,
                    dad  = c(0, 0, 3e9 + 1, 3e9 + 1),
                    mum  = c(0, 0, 3e9 + 2, 3e9 + 2))
  expect_error(build_pedigree(1:3, data = big),
               'Codes out of range: 3000000001, 3000000002, 3000000003, 3000000004$')

  ## The parent columns are checked too
  expect_error(build_pedigree(1:3, data = data.frame(self = 1:2, dad = c(0, 3e9),
                                                     mum = 0)),
               'Codes out of range: 3000000000$')
})


test_that('build_pedigree() stores the recoding map as an integer vector indexed by code (#67)', {

  ## The map is public API: map[x] gives the new codes of original codes x,
  ## and match(y, map) gives back the original codes. So it is an integer
  ## vector as long as the largest code, with NA for codes not in the
  ## pedigree, and no names or other attributes.
  map_of <- function(self, dad, mum)
    attr(suppressWarnings(
      build_pedigree(1:3, data = data.frame(self = self, dad = dad, mum = mum))),
      'map')

  ## Animal 1 precedes its parents 2 and 3 (#49)
  expect_identical(map_of(1:3, c(2L, 0L, 0L), c(3L, 0L, 0L)), c(3L, 1L, 2L))
  ## Each new code is also another animal's original code
  expect_identical(map_of(1:6, c(5L, 5L, 5L, 0L, 0L, 0L), c(6L, 6L, 4L, 0L, 0L, 0L)),
                   c(3L, 4L, 6L, 5L, 1L, 2L))
  ## Consecutive codes starting at 2 (#50)
  expect_identical(map_of(2:4, c(0L, 0L, 2L), c(0L, 0L, 3L)), c(NA, 1L, 2L, 3L))
  ## Gaps, and codes not starting at 1, stored as doubles
  expect_identical(map_of(c(10, 20, 30), c(0, 0, 10), c(0, 0, 20)),
                   replace(rep(NA_integer_, 30), c(10, 20, 30), 1:3))
  ## Offspring coded below its parents, with gaps
  expect_identical(map_of(c(3, 7, 9), c(7, 0, 0), c(9, 0, 0)),
                   c(NA, NA, 3L, NA, NA, NA, 1L, NA, 2L))
  ## Founders added from the parent columns, with an NA parent
  expect_identical(map_of(c(4, 6), c(2, 4), c(1, NA)), c(1L, 2L, NA, 3L, NA, 4L))

  ## The shuffled m4 pedigree above: one slot per code up to the largest,
  ## NA where no animal has that code
  map_fix <- attr(ped_fix, 'map')
  codes <- sort(unique(c(ped_shuffled)))
  codes <- codes[!is.na(codes) & codes > 0]
  expect_identical(typeof(map_fix), 'integer')
  expect_null(attributes(map_fix))
  expect_identical(length(map_fix), as.integer(max(codes)))
  expect_identical(which(!is.na(map_fix)), as.integer(codes))
  expect_identical(sort(map_fix[codes]), seq_along(codes))

  ## A pedigree that needs no recoding has no map
  expect_null(map_of(1:4, c(0, 0, 1, 1), c(0, 0, 2, 2)))
})
