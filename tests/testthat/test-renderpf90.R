
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


test_that('trait masks follow absolute data offsets with and without weights', {
  tr <- c('first', 'middle', 'last')
  base <- effect_group(list(diagonal(factor(rep(1:3, each = 2)))),
                        cov.ini = diag(3), ntraits = 3)
  for (active in list(c(FALSE, TRUE, TRUE), c(TRUE, FALSE, TRUE),
                       c(TRUE, TRUE, FALSE), c(FALSE, TRUE, FALSE))) {
    names(active) <- tr
    group <- effect_group(base$effects, diag(3), 3, trait.active = active)
    for (weighted in c(FALSE, TRUE)) {
      ef <- list(Intercept = fixed(rep(1, 6)), rep = group)
      out <- renderpf90.breedr_modelframe(ef, 3, weighted)
      expected <- rep(5 + weighted, 3)
      expected[!active] <- 0
      expect_identical(out$rep$pos, paste(expected, collapse = ' '))
      expect_identical(out$rep$nest, '')
      expect_identical(out$Intercept$pos,
                       paste(rep(4 + weighted, 3), collapse = ' '))
      expect_identical(out$rep$data, renderpf90(base)$data)
      expect_identical(out$rep$levels, renderpf90(base)$levels)
      expect_identical(dim(out$rep$var), c(3L, 3L))
      expect_true(all(out$rep$var[!active, , drop = FALSE] == 0))
      expect_true(all(out$rep$var[, !active, drop = FALSE] == 0))
    }
  }
})


test_that('all nested virtual positions share the group mask', {
  incidence <- matrix(c(.7,.3,0, 0,.7,.3, .3,0,.7), 3, byrow = TRUE)
  members <- list(a = generic(incidence, covariance = diag(3)),
                   b = generic(incidence, covariance = diag(3)))
  active <- c(first = TRUE, middle = FALSE, last = TRUE)
  group <- effect_group(members, cov.ini = diag(6), ntraits = 3,
                         trait.active = active)
  ef <- list(Intercept = fixed(rep(1, 3)), coupled = group)
  out <- renderpf90.breedr_modelframe(ef, 3, TRUE)$coupled
  pos <- do.call(rbind, strsplit(out$pos, ' '))
  nest <- do.call(rbind, strsplit(out$nest, ' '))
  expect_true(all(pos[, 2] == '0'))
  expect_true(all(nest[, 2] == '0'))
  expect_equal(pos[, 1], pos[, 3])
  expect_equal(nest[, 1], nest[, 3])
  expect_true(all(as.numeric(nest[, 1]) > as.numeric(pos[, 1])))
  expect_identical(dim(out$var), c(6L, 6L))
  expect_true(all(out$var[c(2,5), ] == 0))
  expect_true(all(out$var[, c(2,5)] == 0))
  expect_identical(pf90_effect_layout(ef, 3)$name, c('Intercept', 'a', 'b'))
})


test_that('trait omission preserves structured effect encodings', {
  coord <- as.matrix(expand.grid(x = 1:4, y = 1:4))
  block <- factor(rep(1:4, 4))
  ped <- build_pedigree(1:3, data = data.frame(self = 1:16, sire = 0, dam = 0))
  models <- list(
    diagonal = diagonal(block),
    blocks = breedr_blocks(coord, id = block),
    ar = breedr_ar(coord, rho = c(.3, .4)),
    splines = breedr_splines(coord, n.knots = c(2, 2)),
    animal = additive_genetic_animal(pedigree = ped, idx = 1:16),
    competition = additive_genetic_competition(pedigree = ped,
                    coordinates = coord, id = 1:16, decay = 1, autofill = TRUE),
    pec = permanent_environmental_competition(coordinates = coord,
                                               decay = 1, autofill = TRUE))
  for (nm in names(models)) {
    group <- effect_group(list(models[[nm]]), diag(2), 2)
    masked <- effect_group(list(models[[nm]]), diag(2), 2,
                            trait.active = c(y1 = FALSE, y2 = TRUE))
    reference <- renderpf90.breedr_modelframe(list(group = group), 2, FALSE)$group
    ans <- renderpf90.breedr_modelframe(list(group = masked), 2, FALSE)$group
    expect_equal(ans[c('levels','type','model','file_name','file','data')],
                 reference[c('levels','type','model','file_name','file','data')],
                 label = nm)
    expect_true(all(vapply(strsplit(ans$pos, ' '), `[`, '', 1) == '0'),
                 label = nm)
    nonempty <- nzchar(ans$nest)
    expect_true(all(vapply(strsplit(ans$nest[nonempty], ' '), `[`, '', 1) == '0'),
                 label = nm)
    expect_equal(ans$var, diag(c(0, 1)), label = nm)
  }
})


test_that('parameter groups expose masks only for restricted fits', {
  data <- data.frame(y1 = 1:6, y2 = 6:1, rep = factor(rep(1:3, 2)))
  mc <- call('remlf90', fixed = cbind(y1, y2) ~ 1, random = ~rep,
             data = quote(data))
  mf <- build.mf(mc)
  args <- list(mf, NULL, NULL, NULL, list(rep = diag(2)))
  unrestricted <- do.call(build.effects, args)
  restricted <- do.call(build.effects,
                        c(args, list(traits = list(rep = c(y1 = TRUE, y2 = FALSE)))))
  base <- progsf90(mf, NULL, unrestricted, res.var.ini = diag(2))
  masked <- progsf90(mf, NULL, restricted, res.var.ini = diag(2))
  expect_null(base$parameter$rangroup$rep$trait.active)
  expect_identical(masked$parameter$rangroup$rep$trait.active,
                   c(y1 = TRUE, y2 = FALSE))
  expect_identical(masked$data, base$data)
  expect_identical(masked$files, base$files)
  expect_identical(masked$parameter$rangroup$rep$pos,
                   base$parameter$rangroup$rep$pos)
  expect_identical(masked$parameter$effects$rep, '4 0 3 cross ')
})
