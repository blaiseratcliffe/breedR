
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


test_that("an AR rho grid works when spatial is a variable or forwarded through ... (#5)", {

  ## Each rho of a grid is fitted by re-invoking remlf90(). That used to
  ## rewrite the spatial expression of the matched call, mc$spatial$rho <- rho,
  ## which only works when spatial = list(...) is typed in the call itself: a
  ## variable is a symbol there, and an argument forwarded through ... is ..1,
  ## and neither can be subset. Stop each per-rho fit once its arguments are
  ## checked, and record the rho it was given and its matched call, which is
  ## the call a winning fit stores. No binaries are needed.
  seen <- list()
  calls <- list()
  local_mocked_bindings(
    check_progsf90 = function(...) TRUE,
    build.effects = function(mf, genetic, spatial, ...) {
      seen[[length(seen) + 1]] <<- spatial$rho
      calls[[length(calls) + 1]] <<- get("mcout", envir = parent.frame())
      stop("per-rho fit reached")
    },
    .package = "breedR")

  dat <- expand.grid(x = 1:6, y = 1:5)
  dat$z <- sin(seq_len(nrow(dat)))
  grid <- rbind(c(.2, .5), c(.7, .3))
  rows <- list(c(.2, .5), c(.7, .3))
  sp <- list(model = "AR", coordinates = dat[, c("x", "y")], rho = grid)
  fit_ar <- function(...) remlf90(fixed = z ~ 1, data = dat, ...)
  quietly <- function(x) suppressWarnings(suppressMessages(x))

  ## spatial given as a variable
  expect_error(quietly(remlf90(fixed = z ~ 1, data = dat, spatial = sp)),
               "All rho combinations failed")
  expect_equal(seen, rows)
  ## Each per-rho call carries the evaluated spatial list with its rho set
  expect_identical(lapply(calls, function(cl) cl$spatial),
                   lapply(rows, function(r)
                     list(model = "AR", coordinates = dat[, c("x", "y")],
                          rho = r)))

  ## spatial = list(...) forwarded through a wrapper's ... (the reported case)
  seen <- list()
  expect_error(quietly(fit_ar(spatial = list(model = "AR",
                                             coordinates = dat[, c("x", "y")],
                                             rho = grid))),
               "All rho combinations failed")
  expect_equal(seen, rows)

  ## rho left unset, the default grid, forwarded
  old_ar_eval <- breedR.getOption("ar.eval")
  on.exit(breedR.setOption("ar.eval", old_ar_eval), add = TRUE)
  breedR.setOption("ar.eval", c(-.4, .6))
  seen <- list()
  expect_error(quietly(fit_ar(spatial = sp[c("model", "coordinates")])),
               "All rho combinations failed")
  expect_equal(seen, list(c(-.4, -.4), c(.6, -.4), c(-.4, .6), c(.6, .6)))

  ## The parallel search seeds from its first rho through the same
  ## re-invocation. The mocked seed fit leaves no directory, so this stops at
  ## the seed-directory guard, before makeCluster(): no cluster is started. If
  ## that guard ever moves after makeCluster(), this would start one.
  seen <- list()
  expect_error(quietly(fit_ar(spatial = sp, parallel = 2)),
               "produced no fit directory")
  expect_equal(seen, rows[1])
})


test_that("a literal spatial = list(...) grid builds the same per-rho calls as before (#5)", {

  ## Guard for the fix to #5: when spatial is typed as list(...) in the call,
  ## each per-rho call is still the typed call with its rho element replaced,
  ## so the call a winning fit stores is unchanged. Record the matched call of
  ## each per-rho fit and stop it there. No binaries are needed.
  calls <- list()
  local_mocked_bindings(
    check_progsf90 = function(...) TRUE,
    build.effects = function(mf, genetic, spatial, ...) {
      calls[[length(calls) + 1]] <<- get("mcout", envir = parent.frame())
      stop("per-rho fit reached")
    },
    .package = "breedR")

  dat <- expand.grid(x = 1:6, y = 1:5)
  dat$z <- sin(seq_len(nrow(dat)))
  grid <- rbind(c(.2, .5), c(.7, .3))
  quietly <- function(x) suppressWarnings(suppressMessages(x))

  expect_error(quietly(remlf90(fixed = z ~ 1, data = dat,
                               spatial = list(model = "AR",
                                              coordinates = dat[, c("x", "y")],
                                              rho = grid))),
               "All rho combinations failed")
  expected <- lapply(list(c(.2, .5), c(.7, .3)), function(r)
    bquote(remlf90(fixed = z ~ 1,
                   spatial = list(model = "AR",
                                  coordinates = dat[, c("x", "y")],
                                  rho = .(r)),
                   data = dat)))
  expect_identical(calls, expected)

  ## The parallel seed call, likewise; it stops before any cluster is started
  calls <- list()
  expect_error(quietly(remlf90(fixed = z ~ 1, data = dat,
                               spatial = list(model = "AR",
                                              coordinates = dat[, c("x", "y")],
                                              rho = grid),
                               parallel = 2)),
               "produced no fit directory")
  expect_identical(calls, list(
    bquote(remlf90(fixed = z ~ 1,
                   spatial = list(model = "AR",
                                  coordinates = dat[, c("x", "y")],
                                  rho = .(c(.2, .5))),
                   data = dat, parallel = FALSE))))
})


test_that("select_best_rho() refuses to pick from an all-NA grid", {

  ## A rho can run to completion without error yet still carry no usable
  ## log-likelihood (e.g. an unparseable REML log). which.max(c(NA, NA)) is
  ## integer(0), and ans.rho[[integer(0)]] used to throw "attempt to select
  ## less than one element in get1index" -- an opaque crash several frames
  ## away from the real cause (issue #3).
  expect_error(select_best_rho(c(NA_real_, NA_real_)),
               "usable log-likelihood")
  expect_error(select_best_rho(numeric(0)),
               "usable log-likelihood")

  ## A partially failed grid still picks the best of what succeeded.
  expect_equal(select_best_rho(c(NA, 3.1, -5, 3.9, NA)), 4)
})


test_that("build.AR.rho.grid() honours the ar.eval option", {
  old_ar_eval <- breedR.getOption("ar.eval")
  on.exit(breedR.setOption("ar.eval", old_ar_eval), add = TRUE)
  custom <- c(-.5, -.1, .1, .5)
  breedR.setOption("ar.eval", custom)
  grid <- build.AR.rho.grid(matrix(c(NA_real_, NA_real_), 1, 2))
  expect_equal(sort(unique(grid$rho_r)), sort(custom))
  expect_equal(sort(unique(grid$rho_c)), sort(custom))
})
