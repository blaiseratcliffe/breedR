context("integration helper-testdata guard")

test_that("integration helper-testdata.R does not refit when fixtures already exist", {
  helper_path <- testthat::test_path("..", "integration", "helper-testdata.R")
  skip_if_not(file.exists(helper_path), "tests/integration/helper-testdata.R not found")

  testdata_dir <- system.file("testdata", package = "breedR")
  fixture_names <- c("fixonly", "blk", "ar", "spl", "ped_ar")
  fixture_files <- file.path(testdata_dir, paste0("res_", fixture_names, ".rds"))
  skip_if_not(all(file.exists(fixture_files)),
              "not all integration fixtures are present in inst/testdata")

  ## Snapshot options the helper might touch, and restore them regardless of
  ## how the sys.source() call below exits.
  old_options <- options()
  on.exit(options(old_options), add = TRUE)

  ## Source the helper into an isolated environment where remlf90() is
  ## stubbed to fail loudly. If the helper still tries to fit models when the
  ## fixtures already exist on disk, this stub is what will be called.
  guard_env <- new.env(parent = globalenv())
  guard_env$remlf90 <- function(...) {
    stop("remlf90 should not be called when fixtures already exist")
  }

  expect_no_error(
    sys.source(helper_path, envir = guard_env)
  )
})
