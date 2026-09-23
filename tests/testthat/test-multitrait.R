context("Basic Multitrait models")

test_that("restricted methods distinguish absent effects from unknown estimates", {
  d <- data.frame(y1 = 1:6, y2 = c(3, 2, 4, 6, 5, 7),
                  group = factor(c('b', 'a', 'c', 'b', 'a', 'c')))
  data <- d
  mf <- build.mf(call('remlf90', fixed = cbind(y1, y2) ~ 1,
                      random = ~ group, data = quote(d)))
  effects <- build.effects(mf, NULL, NULL, NULL, list(group = diag(2)),
                            traits = list(group = c(y1 = TRUE, y2 = FALSE)))
  coefficients <- function(value, labels)
    data.frame(value = value, s.e. = ifelse(is.na(value), NA_real_, 1),
                row.names = labels)
  fit <- structure(list(mf = mf, effects = effects, reml = list(method = 'ai'),
                         components = list(pedigree = FALSE),
                         fixed = list(Intercept = list(
                           y1 = coefficients(10, '1'), y2 = coefficients(20, '1'))),
                         ranef = list(group = list(
                           y1 = coefficients(c(-1, 0, 1), levels(d$group)),
                           y2 = coefficients(rep(NA_real_, 3), levels(d$group))))),
                    class = 'remlf90')
  ## Coefficient presence comes from the model, never inferred from NA values.
  expect_equal(unname(fitted(fit)), cbind(10 + c(0, -1, 1, 0, -1, 1), rep(20, 6)))
  expect_identical(dimnames(fitted(fit)), dimnames(as.matrix(model.response(mf))))
  expect_identical(attr(ranef(fit)$group, 'trait.active'), c(y1 = TRUE, y2 = FALSE))
  expect_true(all(is.na(attr(ranef(fit)$group, 'se')[, 'y2'])))
  expect_error(plot(fit), 'Select one trait')
  expect_error(plot(ranef(fit)), 'Select one trait')
  expect_error(vcov(fit), 'does not support trait-restricted')
  fit$ranef$group$y1$value[1] <- NA_real_
  expect_true(anyNA(fitted(fit)[, 'y1']))
  expect_false(anyNA(fitted(fit)[, 'y2']))
})

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
  ## Both fixtures have recoded pedigrees (the douglas codes run to 9764), so
  ## the rows are named by the ids the pedigree was given, which the map takes
  ## back from the internal codes. They used to be the internal codes, which
  ## are other animals' ids (#63).
  ped <- get_pedigree(res)
  map <- attr(ped, 'map')
  expect_false(is.null(map))
  lab <- as.character(match(as.integer(ped@label), map))
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



