
#### Context: breedr_ar() ####
context("AR infrastructure")

# Some coordinates, with a whole row and some points missing
# with non-integer and non-positive value
test.pos <- list(x = c(-2:3), y = 5:8/2)
test.coord <- as.matrix(expand.grid(test.pos$x[-2], test.pos$y))[-8, ]
  # plot(test.coord, pch = 19)
test.rho <- c(.6, .8)
reslst <- list(breedr_ar(test.coord, test.rho, TRUE),
               breedr_ar(test.coord, test.rho, FALSE))

check_build.ar.model <- function(x) {
  inc.mat <- model.matrix(x)
  cov.mat <- get_structure(x)
  eff.size <- ifelse(attr(x, 'grid')$autofill,
                     prod(sapply(test.pos, length)),
                     nrow(test.coord)+1)
  
  test_that("breedr_ar() returns a list with the right elements", {
    expect_is(x, c('ar', 'spatial', 'random', 'breedr_effect'))
    expect_equal(length(x), 5)
    expect_equal(names(x$param), 'rho')
    expect_equal(x$param$rho, test.rho)
    # Incidence matrix
    expect_is(inc.mat, 'sparseMatrix') # a permutation Matrix
    expect_equal(nrow(inc.mat), nrow(test.coord))
    # Covariance matrix
    expect_is(cov.mat, 'sparseMatrix') 
    expect_equal(ncol(cov.mat), eff.size)
    expect_equal(ncol(inc.mat), nrow(cov.mat))
  })
}

for (x in reslst) check_build.ar.model(x)


#### Context: AR rho grid search ####
context("AR rho grid search")

## The grid fits one model per rho and picks the most likely. Argument
## combinations that leave it with nothing to rank are refused here, before any
## fitting starts -- these need no binaries for that reason.

test_that("an AR rho grid is refused for a non-local fit", {

  ## Hide the backend, as CI installs without it, so that the claim above is
  ## tested wherever this runs and not only on machines that lack it.
  no_bin <- tempfile("no_progsf90_")
  dir.create(no_bin)
  old_bin <- breedR.getOption("breedR.bin")
  breedR.setOption("breedR.bin", no_bin)
  on.exit({
    breedR.setOption("breedR.bin", old_bin)
    unlink(no_bin, recursive = TRUE)
  }, add = TRUE)
  expect_false(check_progsf90(quiet = TRUE))

  grid_fit <- function(bin, ...)
    remlf90(fixed = phe_X ~ gg, data = globulus,
            spatial = list(model = 'AR',
                           coord = globulus[, c('x', 'y')],
                           rho = rbind(c(.8, .8), c(.9, .9))),
            breedR.bin = bin, ...)

  ## A submitted fit returns a job id rather than a likelihood, so the grid has
  ## no way to rank its rhos. It used to get as far as building a binary path
  ## out of breedR.bin -- literally "submit/blupf90+" -- and report that every
  ## rho had failed.
  msg <- tryCatch(grid_fit('submit'), error = conditionMessage)
  expect_match(msg, "requires a local fit")
  expect_match(msg, "submit", fixed = TRUE)
  expect_false(grepl("All rho combinations failed", msg, fixed = TRUE))

  ## 'remote' is the same story with the results fetched back afterwards.
  expect_error(grid_fit('remote'), "requires a local fit")

  ## The name is matched case-insensitively, as it is everywhere else it is
  ## tested (see check_progress_args()).
  expect_error(grid_fit('Submit'), "requires a local fit")

  ## The log and resume refusal beside it is an argument check too.
  expect_error(grid_fit(no_bin,
                        progress_file = file.path(tempdir(), "grid.log")),
               "rho grid")

  ## A local grid still needs the backend, and must be told so before the
  ## per-rho fits start: each runs inside a tryCatch, which would bury the
  ## cause under "All rho combinations failed".
  expect_error(grid_fit(no_bin), "Binary dependencies missing")
})
