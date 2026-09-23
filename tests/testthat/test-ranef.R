context("Names of the random effects and their standard errors")

## parse_results() keeps the solutions file's row numbers as the names of each
## effect's values and s.e. ranef() relabels the levels (animals, or the rows
## of a generic structure matrix), and must relabel the s.e. alike: otherwise
## a lookup by name returns another level's s.e. (issue #60).

test_that("the s.e. of the genetic effect is named by animal", {

  res <- load_res("ped_ar")
  expect_null(attr(get_pedigree(res), 'map'))   # not recoded: labels are the ids
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


test_that("genetic values are named by the ids the pedigree was given, when it is recoded (#63)", {

  ## Parents 3 and 4 are coded above their offspring 1 and 2, so
  ## build_pedigree() recodes 3, 4, 1, 2 as 1, 2, 3, 4. Every internal code
  ## is another animal's id: named by them, animal 1's value was found
  ## under '3', and animal 3's under '1'.
  ped <- suppressWarnings(build_pedigree(1:3, data = data.frame(
    self = c(3, 4, 1, 2), dad = c(0, 0, 3, 3), mum = c(0, 0, 4, 4))))
  expect_identical(attr(ped, 'map'), c(3L, 4L, 1L, 2L))

  ## The solutions list the animals by internal code: row k is the animal
  ## recoded as k
  id  <- c(2, 1, 2)
  res <- structure(
    list(components = list(pedigree = TRUE),
         effects = list(genetic = effect_group(
           list(direct = additive_genetic_animal(ped, id)),
           cov.ini = 1, ntraits = 1)),
         ranef = list(genetic = list(data.frame(
           value = c(0.1, 0.2, 0.3, 0.4), s.e. = c(1, 2, 3, 4),
           row.names = 1:4)))),
    class = 'remlf90')

  g  <- ranef(res)$genetic
  se <- attr(g, 'se')
  expect_identical(names(g), c('3', '4', '1', '2'))
  expect_identical(names(se), names(g))

  ## renamed only: the same numbers in the same order
  expect_identical(as.numeric(g),  c(0.1, 0.2, 0.3, 0.4))
  expect_identical(as.numeric(se), c(1, 2, 3, 4))

  ## each id finds its own animal
  expect_identical(g[['1']], 0.3)
  expect_identical(se[['3']], 1)

  ## and the column each record points to is named by the record's id
  Z <- as.matrix(res$effects$genetic$effects$direct$incidence.matrix)
  expect_identical(names(g)[apply(Z, 1, which.max)], as.character(id))
})

test_that("pedigree_labels() translates recoded codes back, and leaves others alone (#63)", {
  p0 <- build_pedigree(1:3, data = data.frame(self = 1:4, dad = c(0, 0, 1, 1),
                                              mum = c(0, 0, 2, 2)))
  expect_null(attr(p0, 'map'))
  expect_identical(pedigree_labels(p0), p0@label)

  ## gaps, and offspring coded below a parent
  p1 <- suppressWarnings(build_pedigree(1:3, data = data.frame(
    self = c(30, 7, 12, 5), dad = c(0, 0, 30, 12), mum = c(0, 0, 7, 7))))
  map <- attr(p1, 'map')
  expect_identical(pedigree_labels(p1),
                   as.character(match(as.integer(p1@label), map)))
  expect_identical(map[as.integer(pedigree_labels(p1))], seq_along(p1@label))
})
