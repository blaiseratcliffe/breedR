context("Basic Multitrait models")

## Bivariate model with a random effect of block and a genetic effect.
res_mt <- readRDS(file.path(testdata, "res_mt.rds"))
ntraits <- ncol(as.matrix(model.response(res_mt$mf)))

test_that("Extract random effects", {
  
  ## ranef() recovers a list of all existing random effects
  expect_identical(length(ranef(res_mt)), length(res_mt$ranef))
  
  ## each element is a matrix with as many columns as traits
  expect_identical(lapply(ranef(res_mt), ncol), lapply(res_mt$ranef, length))
})


## The rows of a multi-trait genetic effect, and of its s.e., are its animals.
## They used to keep the row numbers of the solutions file, so a lookup by
## animal returned another animal's breeding values (issue #60).
expect_genetic_rows_by_animal <- function(res) {
  lab <- as.character(get_pedigree(res)@label)
  rr <- ranef(res)
  gen.idx <- grep('genetic', names(rr))
  expect_true(length(gen.idx) > 0)
  for (k in gen.idx) {
    g <- rr[[k]]
    se <- attr(g, 'se')
    expect_identical(rownames(g), lab)
    expect_identical(rownames(se), lab)
    expect_null(names(g))

    ## the solutions list the animals in code order, one data frame per trait
    i <- 105
    expect_equal(unname(g[lab[i], ]),
                 vapply(res$ranef[[k]], function(d) d$value[i], 0,
                        USE.NAMES = FALSE))
    expect_equal(unname(se[lab[i], ]),
                 vapply(res$ranef[[k]], function(d) d$s.e.[i], 0,
                        USE.NAMES = FALSE))
  }
}

test_that("Genetic values and their s.e. are labelled by animal", {
  expect_genetic_rows_by_animal(res_mt)
})



context("Multitrait-competition models")

## Bivariate model with a random effect of block and a genetic effect.
res_mtcp <- readRDS(file.path(testdata, "res_mtcp.rds"))
ntraits <- ncol(as.matrix(model.response(res_mtcp$mf)))

test_that("Extract random effects", {
  
  ## ranef() recovers a list of all existing random effects
  expect_identical(length(ranef(res_mtcp)), length(res_mtcp$ranef))
  
  ## each element is a matrix with as many columns as traits
  expect_identical(lapply(ranef(res_mtcp), ncol), lapply(res_mtcp$ranef, length))
})

test_that("Direct and competition values and their s.e. are labelled by animal", {
  expect_genetic_rows_by_animal(res_mtcp)
})



