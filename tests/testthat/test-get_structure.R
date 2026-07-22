### Test the functions for getting the structure matrices ###


#### Context: Extracting structure matrices ####
context("Extracting structure matrices")

## Extracting structure matrices from simple breedR effects

coord <- expand.grid(list(x = seq(1, 100, length = 51),
                          y = seq(1001, 1100, length = 35)),
                     KEEP.OUT.ATTRS = FALSE)

# A splines object
spl <- breedr_splines(coord)
spl.str <- get_structure(spl)

# A generic object (with same structure, but inverted)
inv_spl.str <- solve(spl.str)
gen <- generic(incidence = model.matrix(spl),
               precision = inv_spl.str)
gen.str <- get_structure(gen)

test_that('get_structure() extracts a Matrix', {
  expect_is(spl.str, 'Matrix')
  expect_is(gen.str, 'Matrix')
})


test_that('get_structure() recovers the right structure type', {
  expect_identical(attr(spl.str, 'type'), 'covariance')
  expect_identical(attr(gen.str, 'type'), 'precision')
})


## Extracting structure matrices from groups of effects

eg <- effect_group(list(spl, gen), cov.ini = diag(1,2,2), ntraits = 1)
eg.str <- get_structure(eg)

test_that('get_structure() recovers the common structure in Matrix format', {
  expect_is(eg.str, 'Matrix')
  # Compare dense content, not the S4 object: solve() above populates spl.str's
  # cached @factors slot, which is irrelevant to the structure but breaks a
  # strict object comparison.
  expect_equal(as.matrix(eg.str), as.matrix(spl.str))
})


## A group of 3+ effects sharing the same structure type must not crash.
## Regression test: str.list[[-1]] is an invalid subscript for length >= 3.
test_that('get_structure() handles groups of 3+ same-type effects', {
  spl.b <- breedr_splines(coord)
  spl.c <- breedr_splines(coord)
  eg3 <- effect_group(list(spl, spl.b, spl.c), cov.ini = diag(1, 3, 3),
                      ntraits = 1)
  expect_error(eg3.str <- get_structure(eg3), NA)  # no error
  expect_is(eg3.str, 'Matrix')
  expect_equal(as.matrix(eg3.str), as.matrix(spl.str))
})


## Extracting structure matrices from breedR objects

test_that('get_structure() retrieves an empty list from a model fit without random effects', {

  res <- load_res("fixonly")
  breedr.str <- get_structure(res)
  
  expect_is(breedr.str, 'list')
  expect_equal(breedr.str, list(), check.attributes = FALSE)
})


test_that('get_structure() retrieves a list of structure matrices from a model fit', {

  res <- load_res("ar")
  breedr.str <- get_structure(res)
  
  expect_is(breedr.str, 'list')
  for (i in seq_along(breedr.str)) {
    expect_is(breedr.str[[i]], "Matrix")
  }
  
})
