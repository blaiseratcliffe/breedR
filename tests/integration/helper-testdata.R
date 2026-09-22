# Generate (or update) testdata used by unit tests

## target directory
testdata <- system.file("testdata", package = "breedR")

## Only (re)generate the fixtures when at least one is missing. reml$output
## carries CPU timing lines (see progsf90.R), so refitting on every run makes
## these committed fixtures differ on every run even though nothing about the
## model changed. Skip the fit (and the save) entirely once all 5 exist.
target_files <- file.path(testdata,
                           paste0("res_", c("fixonly", "blk", "ar", "spl", "ped_ar"), ".rds"))
if (!all(file.exists(target_files))) {
  message("Generating test data ...")

  #### fitted models ####

  res  <- suppressMessages(
    list(
      ## A simple model with one fixed effect
      fixonly = remlf90(
        fixed  = phe_X ~ gg,
        data = globulus
      ),
      ## A spatial blocks model
      blk = remlf90(
        fixed  = phe_X ~ 1,
        spatial = list(
          model = 'blocks',
          coord = globulus[, c('x','y')],
          id = 'bl'
        ),
        data = globulus
      ),
      ## An spatial autoregressive model with one random effect
      ar = remlf90(
        fixed  = phe_X ~ 1,
        random = ~ gg,
        spatial = list(
          model = 'AR',
          coord = globulus[, c('x','y')],
          rho = c(.85, .8)
        ),
        data = globulus
      ),
      ## An spatial splines (2x2 knots) model fitted with EM
      spl = remlf90(
        fixed  = phe_X ~ 1,
        spatial = list(
          model = 'splines',
          coord = globulus[, c('x','y')],
          n.knots = c(2, 2)
        ),
        data = globulus,
        method = 'em'
      ),
      ## A genetic-AR model with a fixed effect
      ped_ar = remlf90(
        fixed  = phe_X ~ gg,
        genetic = list(
          model = 'add_animal',
          pedigree = globulus[,1:3],
          id = 'self'
        ),
        spatial = list(
          model = 'AR',
          coord = globulus[, c('x','y')],
          rho = c(.85, .8)
        ),
        data = globulus
      )
    )
  )


  ## Drop the working directory before saving. It is an absolute path under the
  ## tempdir() of whoever last ran this, so leaving it in makes every rerun a
  ## diff against the committed fixture and gives the next reader a path that
  ## went away with a session they never had. Nothing in the unit tests reads it.
  for (idx in seq_along(res)){
    fn <- paste0("res_", names(res)[idx], ".rds")
    fit <- res[[idx]]
    fit$reml$dir <- NULL
    saveRDS(fit, file = file.path(testdata, fn))
  }
}

