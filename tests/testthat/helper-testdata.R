## test data directory
testdata <- system.file("testdata", package = "breedR")

# Quick convenience function to fit a test model on globulus data
breedR.result <- function(...) {
  dat <- breedR::globulus
  suppressMessages(
    remlf90(fixed  = phe_X ~ gg,
            genetic = list(model = 'add_animal',
                           pedigree = dat[, 1:3],
                           id = 'self'),
            spatial = list(model = 'AR',
                           coord = dat[, c('x', 'y')],
                           rho = c(.85, .8)),
            data = dat,
            ...)
  )
}

load_res <- function(key, dir = testdata) {
  fn <- paste0("res_", key, ".rds")
  readRDS(file.path(dir, fn))
}


## Generate sample REML output files

## 10-variate model.
## large residual covariance matrix (10x10) which wraps lines
if (!file.exists(tf <- file.path(testdata, "airemlf90_log_4.txt"))) {
  
  tdat <- data.frame(replicate(n = 10, rnorm(1e2)))
  res <- remlf90(
    cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) ~ 1,
    data = tdat,
    method = "ai"
  )
  writeLines(res$reml$output, tf)
}


## Generate fitted multi-trait model
if (!file.exists(tf <- file.path(testdata, "res_mt.rds"))) {
  dat <-
    droplevels(
      douglas[douglas$site == "s3",
              names(douglas)[!grepl("H0[^4]|AN|BR|site",
                                    names(douglas))]]
    )
  res_mt <-
    remlf90(
      fixed = cbind(H04, C13) ~ orig,
      random = ~ block,
      genetic = list(
        model = 'add_animal',
        pedigree = dat[, 1:3],
        id = 'self'),
      data = dat
    )

  saveRDS(res_mt, file = file.path(testdata, "res_mt.rds"))
}

## Generate fitted multitrait-competition model
if (!file.exists(tf <- file.path(testdata, "res_mtcp.rds"))) {
  dat <-
    droplevels(
      douglas[douglas$site == "s3",
              names(douglas)[!grepl("H0[^4]|AN|BR|site",
                                    names(douglas))]]
    )
  res_mtcp <-
    remlf90(
      fixed = cbind(H04, C13) ~ orig,
      genetic = list(
        model = 'competition',
        pedigree = dat[, 1:3],
        id = 'self',
        coord = dat[, c("x", "y")],
        pec = FALSE
      ),
      data = dat,
      method = "em"
    )
  
  saveRDS(res_mtcp, file = file.path(testdata, "res_mtcp.rds"))
}

