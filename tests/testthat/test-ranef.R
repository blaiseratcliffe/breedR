context("Names of the random effects and their standard errors")

## parse_results() keeps the solutions file's row numbers as the names of each
## effect's values and s.e. ranef() relabels the levels (animals, or the rows
## of a generic structure matrix), and must relabel the s.e. alike: otherwise
## a lookup by name returns another level's s.e. (issue #60).

test_that("the s.e. of the genetic effect is named by animal", {

  res <- load_res("ped_ar")
  g <- ranef(res)$genetic
  se <- attr(g, 'se')
  lab <- as.character(get_pedigree(res)@label)

  expect_identical(names(g), lab)
  expect_identical(names(se), lab)

  ## the solutions list the animals in code order, so animal 1000's s.e. is
  ## the 1000th of the effect. The name '1000' used to belong to animal 986.
  expect_identical(lab[1000], '1000')
  expect_equal(se[['1000']], res$ranef$genetic[[1]]$s.e.[1000])
})


test_that("the s.e. of a generic effect is named like its values", {

  ## A fit reduced to what ranef() reads: a generic effect whose structure
  ## matrix names its levels, and its parsed solutions, which are named by
  ## their row in the solutions file.
  lv <- c('a', 'b', 'c')
  cov.mat <- diag(3)
  dimnames(cov.mat) <- list(lv, lv)
  res <- structure(
    list(components = list(pedigree = FALSE),
         effects = list(
           g = effect_group(list(generic(incidence = diag(3),
                                         covariance = cov.mat)),
                            cov.ini = 1, ntraits = 1)),
         ranef = list(
           g = list(data.frame(value = c(0.1, 0.2, 0.3),
                               s.e. = c(1, 2, 3),
                               row.names = 5:7)))),
    class = 'remlf90')

  g <- ranef(res)$g
  expect_identical(names(g), lv)
  expect_identical(names(attr(g, 'se')), lv)
  expect_equal(attr(g, 'se')[['b']], 2)
})
