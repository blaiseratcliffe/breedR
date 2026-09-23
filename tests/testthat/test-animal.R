
data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

#### Context: Building additive-genetic models ####
context("Building additive-genetic models")

test_that("Correctly builds structure of additive_genetic_animal component", {
  dat <- data.frame(id = 1:4,
                    sire = c(11, 11, 2, 3),
                    dam  = c(12, NA, 1, 12))
  ## Recoded pedigree with further unobserved descendants
  ## (with will come later in the codification)
  ped <- suppressWarnings(build_pedigree(1:3, data = rbind(dat, c(5, 1, 2))))
  aga <- try(breedR:::additive_genetic_animal(ped, dat$id))
  
  expect_true(!inherits(aga, "try-error"))
  expect_is(aga, 
            c("additive_genetic_animal", 
              "additive_genetic",
              "genetic",
              "random", 
              "breedr_effect"))
  expect_named(aga, 
               c("incidence.matrix", 
                 "structure.matrix",
                 "structure.type", 
                 "pedigree"))
  expect_identical(dim(aga$incidence.matrix),
                   c(nrow(dat), nrow(as.data.frame(ped))))
  expect_identical(dim(aga$structure.matrix),
                   rep(nrow(as.data.frame(ped)), 2))
  expect_identical(aga$structure.type, 'covariance')
  expect_identical(aga$pedigree, ped)
})

test_that("additive_genetic_animal() gives one incidence row per id, or refuses (#62)", {
  ped <- suppressWarnings(build_pedigree(1:3, data = data.frame(
    self = c(10, 20, 30, 40), dad = c(0, 0, 10, 10), mum = c(0, 0, 20, 20))))
  map <- attr(ped, 'map')
  expect_false(is.null(map))

  ## Valid ids, one repeated: one row each, in the column of its own animal
  inc  <- additive_genetic_animal(ped, c(40, 30, 40))$incidence.matrix
  orig <- match(as.integer(ped@label), map)
  expect_identical(nrow(inc), 3L)
  expect_identical(orig[apply(as.matrix(inc), 1, which.max)], c(40L, 30L, 40L))

  ## An id of 0 used to be dropped: 2 rows for 3 records, and no error here
  expect_error(additive_genetic_animal(ped, c(30, 40, 0)),
               'not represented in the pedigree:\n 0', fixed = TRUE)
  expect_error(additive_genetic_animal(ped, c(30, 40.5)),
               'not represented in the pedigree:\n 40.5', fixed = TRUE)
})
