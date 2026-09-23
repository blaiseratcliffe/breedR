

#### Context: competition infrastructure ####
context("competition infrastructure")

## Minimal dataset
dat <- data.frame(id   = 1:6,
                  sire = c(11, 11, 2, 3, 11, 3),
                  dam  = c(12, NA, 1, 12, 12, 1),
                  x    = c(1,2,-1,0,0,1),
                  y    = c(-1,0,0,1,-1,1))
## Corresponding pedigree with additional offspring
ped <- suppressWarnings(build_pedigree(1:3, data = rbind(dat, c(7, 1, 2))))
var.ini.mat <- matrix(c(1, -.5, -.5, 1), 2, 2)


test_that("Valid alternative model specifications pass check_genetic()", {
  ## specify var.ini, but use pec=FALSE (as by default)
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = ped,
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))),
    NA
  )
  
  ## non-recoded pedigree
  idx <- attr(ped, 'map')[dat$id]
  dat[, 1:3] <- as.data.frame(ped)[idx, ]
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = dat[, 1:3],
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))),
    NA
  )
  
})



test_that("Invalid alternative model specifications fail check_genetic()", {

  ## incomplete non-recoded pedigrees
  idx <- attr(ped, 'map')[dat$id]
  dat[, 1:3] <- as.data.frame(ped)[idx, ]
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = dat[-nrow(dat), 1:3],
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))), 
    'The following individuals in id are not represented'
  )
})



test_that("additive_genetic_competition() works as expected", {
  ## Full specification, from minimal input
  comp.spec <- check_genetic(
    model = 'competition',
    pedigree = ped,
    id = dat$id,
    coordinates = dat[, c('x', 'y')],
    pec = TRUE,
    response = rnorm(nrow(dat))
  )
  
  res <- with(
    comp.spec,
    additive_genetic_competition(
      pedigree    = pedigree,
      coordinates = coordinates,
      id          = id,
      decay       = competition_decay,
      autofill    = autofill
    )
  )
  
  expect_is(res, c("additive_genetic_competition", "additive_genetic", "genetic", 
                   "competition", "spatial", "random", "breedr_effect"))
  expect_equal(length(res), 5)
  # Incidence matrix
  inc.mat <- model.matrix(res)
  expect_is(inc.mat, 'sparseMatrix') # a permutation Matrix
  expect_equal(nrow(inc.mat), nrow(dat))
  # Covariance matrix
  cov.mat <- get_structure(res)
  expect_is(cov.mat, 'sparseMatrix') 
  expect_equal(ncol(cov.mat), nrow(as.data.frame(ped)))
  expect_equal(ncol(inc.mat), nrow(cov.mat))
})

test_that("the competition incidence follows the recode map (#65), and an id of 0 is refused (#62)", {
  ## 4 founders (1-4) and 12 trees (5-16) on a regular 4 x 3 grid, coded
  ## 1..16 in order: this pedigree is not recoded.
  ped1  <- data.frame(self = 1:16,
                      dad  = c(0, 0, 0, 0, rep(1:2, 6)),
                      mum  = c(0, 0, 0, 0, rep(3:4, each = 6)))
  trees <- data.frame(id = 5:16, x = rep(1:4, 3), y = rep(1:3, each = 4))

  ## The same trial relabelled: gaps, and most trees coded below their
  ## parents, so build_pedigree() recodes and reorders it.
  new_code <- c(900, 850, 870, 999, 20, 5, 710, 33, 150, 1, 64, 400, 12, 300, 77, 2)
  relabel  <- function(x) c(0, new_code)[x + 1]    # 0 (unknown parent) stays 0
  ped2 <- as.data.frame(lapply(ped1, relabel))

  spec <- function(pedigree, id)
    suppressWarnings(check_genetic(model = 'competition', pedigree = pedigree,
                                   id = id, coordinates = trees[, c('x', 'y')],
                                   var.ini = var.ini.mat,
                                   response = seq_len(nrow(trees))))
  comp_inc <- function(s)
    as.matrix(additive_genetic_competition(s$pedigree, s$coordinates, s$id,
                                           s$competition_decay,
                                           s$autofill)$incidence.matrix)
  s1 <- spec(ped1, trees$id)
  s2 <- spec(ped2, relabel(trees$id))
  expect_null(attr(s1$pedigree, 'map'))
  expect_false(is.null(map2 <- attr(s2$pedigree, 'map')))

  ## Original ids of the columns of the recoded fit, in ped1's coding
  col_id <- match(match(as.integer(s2$pedigree@label), map2), new_code)
  expect_setequal(col_id, 1:16)

  ## Each record, and each of its neighbours, sits on the same animal in both
  inc1 <- comp_inc(s1)
  inc2 <- comp_inc(s2)
  expect_identical(dim(inc2), c(12L, 16L))
  expect_identical(inc2[, order(col_id)], inc1)
  dir1 <- as.matrix(additive_genetic_animal(s1$pedigree, s1$id)$incidence.matrix)
  dir2 <- as.matrix(additive_genetic_animal(s2$pedigree, s2$id)$incidence.matrix)
  expect_identical(dir2[, order(col_id)], dir1)

  ## An id of 0 used to pass check_genetic() and then shrink the lookup, which
  ## failed here with "number of items to replace ..."
  bad <- replace(relabel(trees$id), 12, 0)
  msg <- 'not represented in the pedigree:\n 0'
  expect_error(spec(ped2, bad), msg, fixed = TRUE)
  expect_error(additive_genetic_competition(s2$pedigree, s2$coordinates, bad,
                                            s2$competition_decay, s2$autofill),
               msg, fixed = TRUE)
})



test_that("neighbours.at.list() accepts a list of matrices (R 4.0 class regression, #27)", {
  m1 <- matrix(1:4, 2)
  m2 <- matrix(5:8, 2)
  expect_error(neighbours.at(list(a = m1, b = m2), "N"), NA)
  res <- neighbours.at(list(a = m1, b = m2), "N")
  expect_equal(res$a, neighbours.at(m1, "N"))
  expect_equal(res$b, neighbours.at(m2, "N"))

  ## Multiple directions (the sapply branch of neighbours.at.matrix())
  dirs <- c("N", "S", "E", "W")
  res.dirs <- neighbours.at(list(a = m1, b = m2), dirs)
  expect_equal(res.dirs$a, neighbours.at(m1, dirs))
  expect_equal(res.dirs$b, neighbours.at(m2, dirs))

  ## The guard still rejects a non-matrix element
  expect_error(neighbours.at(list(a = m1, b = 1:4), "N"))
})
