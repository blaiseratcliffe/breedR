
### Test the rendering to pf90 format ###

context("Render PF90")

test_that('renderpf90.matrix() renders different matrix types', {

  mat.list <- list(
    ## one effective column, some rows full-zero
    list(m = diag(c(2, 0, 4)),
         r = matrix(c(2, 0, 4, 1, 0, 3), ncol = 2))
    ,
    ## two effective columns, some rows full-zero
    ## some column fill-in with zero also needed
    list(m = rbind(c(11, 13, 0), 0, c(0, 0, 14)),
         r = rbind(c(11, 13, 1, 2), 0, c(14, 0, 3, 0)))
  )
  
  for (x in mat.list) {
    res <- try(renderpf90.matrix(x$m))
    
    expect_true(!(failed <- inherits(res, 'try-error')))
    
    if (!failed) {
      expect_equal(res, x$r)
    }
  }

})


test_that('renderpf90.breedr_modelframe() renders a single-trait breedr_modelframe correctly', {

  # TODO...  
  # testdat <- transform(
  #   expand.grid(x = 1:4, y = 1:4, KEEP.OUT.ATTRS = FALSE),
  #   z = rnorm(16),
  #   mu = 1)
  # 
  # bc <- call('remlf90', fixed = phe_X~1, data = quote(as.data.frame(m1)))
  # str(build.mf(bc))
  # 
  # 
  # breedrmf <- build.effects(mf = build.mf(bc),
  #                           genetic = NULL,
  #                           spatial = NULL,
  #                           generic = NULL)
  # renderpf90.breedr_modelframe(breedrmf, ntraits = 1)
})


test_that('as.triplet() keeps an implicit diagonal', {

  ## #22. A unit-diagonal Matrix stores nothing in @x -- its diagonal is
  ## recorded in a flag -- so extracting the stored entries returned nothing
  ## and the structure file was written empty. The backend read 0 elements,
  ## produced no solutions, and the fit failed several steps later on a file
  ## that had never been written.
  n <- 4

  ## every spelling of an identity that can reach here
  identities <- list(
    ddi       = Matrix::Diagonal(n),
    as.Matrix = as(as.matrix(diag(n)), 'Matrix'),   # what random() produces
    base      = diag(n),
    dgC       = as(as(as(diag(n), 'dMatrix'), 'generalMatrix'), 'CsparseMatrix'),
    ind       = as(seq_len(n), 'indMatrix')         # no @x slot at all
  )

  for (nm in names(identities)) {
    tri <- as.triplet(identities[[nm]])
    expect_equal(nrow(tri), n, label = paste('rows for', nm))
    expect_equal(ncol(tri), 3L, label = paste('columns for', nm))
    ## identity: the diagonal, all ones
    expect_equal(tri[, 1], tri[, 2], label = paste('diagonal for', nm))
    expect_equal(tri[, 3], rep(1, n), label = paste('values for', nm))
  }
})


test_that('as.triplet() gives the lower triangle of a symmetric matrix', {

  ## Row *order* is deliberately not asserted. Coercing through generalMatrix
  ## to make the diagonal explicit also makes the entries come out
  ## column-major, where the single-triangle classes used to give row-major.
  ## Same set either way, and the backend reads triplets in any order.
  m <- Matrix::Matrix(c(2, 1, 0, 0,
                        1, 2, 1, 0,
                        0, 1, 2, 1,
                        0, 0, 1, 2), nrow = 4, sparse = TRUE)
  m <- as(as(m, 'symmetricMatrix'), 'CsparseMatrix')

  tri <- as.triplet(m)

  ## the lower triangle of the full matrix, and nothing above it
  expect_true(all(tri[, 1] >= tri[, 2]))
  expect_equal(nrow(tri), sum(as.matrix(m)[lower.tri(m, diag = TRUE)] != 0))

  ## every stored entry is the value the matrix actually holds
  expect_equal(tri[, 3], as.matrix(m)[cbind(tri[, 1], tri[, 2])])
})


test_that('renderpf90.generic() writes a structure file for an identity', {

  ## The end of the #22 path that does not need the binaries: whatever
  ## as.triplet() returns is written verbatim by write.progsf90(), so a
  ## zero-row triplet means a zero-byte file with certainty.
  n <- 6
  inc <- as(rep(seq_len(n), 2), 'indMatrix')   # 2n x n incidence

  gm <- generic(incidence = inc, covariance = diag(n))
  ans <- renderpf90(gm)

  expect_equal(nrow(ans$file), n)
  expect_equal(ans$file_name, 'generic')
  ## a covariance, so the backend is asked to invert it
  expect_equal(ans$model, 'user_file_i')
})
