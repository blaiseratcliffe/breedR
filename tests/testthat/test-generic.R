
### Test the generic model ###

context("Generic infrastructure")

test.inc8x4 <- as(rep(1:4, 2), 'indMatrix')   # An 8x4 incidence matrix
test.cov4x4 <- with(L <- Matrix::tril(matrix(sample(16),4)),
                    Matrix::t(L)%*%L)         # a SPD 4x4 covariance matrix

test_that('generic() takes an incidence matrix in all matrix-like formats', {
  
  ## logical matrix (converted to Matrix::lgCMatrix)
  logical.matrix <- as.matrix(test.inc8x4)
  gm1 <- generic(incidence  = logical.matrix, covariance = test.cov4x4)
  expect_identical(gm1$incidence.matrix, as(logical.matrix, 'Matrix'))
  
  ## numeric dataframe (convert to matrix and then to Matrix::dgCMatrix)
  numeric.dataframe <- with(list(logical.matrix), {
    mode(logical.matrix) <- 'double'
    data.frame(logical.matrix)})
  gm2 <- generic(incidence  = numeric.dataframe, covariance = test.cov4x4)
  expect_identical(gm2$incidence.matrix,
                   as(as.matrix(numeric.dataframe), 'Matrix'))
  
  ## index Matrix (not converted at all)
  gm3 <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
  expect_identical(gm3$incidence.matrix, test.inc8x4)
  
})


test_that('generic() takes either covariance or precision matrices', {
  
  gm21 <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
  expect_match(gm21$structure.type, 'covariance')
  
  gm22 <- generic(incidence  = test.inc8x4, precision = test.cov4x4)
  expect_match(gm22$structure.type, 'precision')

})


test_that('generic() stores an identity covariance as a unit-diagonal Matrix', {

  ## Pins the shape that caused #22, so the class the fix has to cope with is
  ## on the record rather than an incidental detail of Matrix's constructor.
  ## `diag(n)` is the obvious covariance to hand to generic(), and as.Matrix()
  ## turns it into a ddiMatrix whose diagonal is implicit -- length(@x) is 0,
  ## and extracting the stored entries used to yield nothing at all.
  gm <- generic(incidence = test.inc8x4, covariance = diag(4))

  expect_s4_class(gm$structure.matrix, 'ddiMatrix')
  expect_identical(gm$structure.matrix@diag, 'U')
  expect_length(gm$structure.matrix@x, 0L)

  ## and the rendering survives it
  expect_equal(nrow(as.triplet(vcov(gm))), 4L)
})

## TODO:
## Test more throughly:
#   - effect_type()
#   - renderpf90()
#   - all.equal.remlf90()
#   - breedr_effect()
#   - model.matrix.effect_group()
#   - model.matrix.breedr_effect()
#   - random()
#   - vcov.random()
