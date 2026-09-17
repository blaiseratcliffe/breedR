## Tests for heterogeneous residual variance helpers (R/hetres.R)

context("Heterogeneous residual variance helpers")

test_that("hetres_options produces class-based options", {
  result <- hetres_options(group_col = 5, n_groups = 10, var_file = "hetres")
  expect_match(result[1], "^hetres_int 5 10$")
})

test_that("hetres_options requires var_file for class-based residuals", {
  expect_error(hetres_options(group_col = 5, n_groups = 10),
               "var_file.*required")
})

test_that("hetres_options produces covariate-based options", {
  result <- hetres_options(covariate_cols = c(8, 9),
                            initial = c(4.0, 0.1, 0.1))
  expect_length(result, 2)
  expect_match(result[1], "^hetres_pos 8 9$")
  expect_match(result[2], "^hetres_pol 4 0\\.1 0\\.1$")
})

test_that("hetres_options includes var_file for class-based", {
  result <- hetres_options(group_col = 5, n_groups = 3,
                            var_file = "my_hetres")
  expect_length(result, 2)
  expect_match(result[2], "^hetres_var my_hetres$")
})

test_that("hetres_options works with single covariate", {
  result <- hetres_options(covariate_cols = 10,
                            initial = c(2.5, 0.01))
  expect_match(result[1], "^hetres_pos 10$")
  expect_match(result[2], "^hetres_pol 2\\.5 0\\.01$")
})

test_that("hetres_options covariate-based without initial values", {
  result <- hetres_options(covariate_cols = 8)
  expect_length(result, 1)
  expect_match(result[1], "^hetres_pos 8$")
})

test_that("hetres_options errors with no arguments", {
  expect_error(hetres_options(), "Specify either")
})

test_that("hetres_options errors with both group and covariate", {
  expect_error(
    hetres_options(group_col = 5, n_groups = 3, covariate_cols = 8),
    "not both"
  )
})

test_that("hetres_options errors when n_groups missing for class-based", {
  expect_error(
    hetres_options(group_col = 5),
    "n_groups.*required"
  )
})


context("Class-based heterogeneous residuals under REML")

## BLUPF90+ 2.73 reads OPTION hetres_int and its variance file under
## 'method VCE', echoes the class variances as "Fixed R", and then fits and
## solves a homoscedastic model: under AI and EM-REML alike, the estimates,
## -2logL and solutions match a fit without the option to every digit (#39).
## Every breedR REML entry point writes 'method VCE', so each one refuses it.
##
## The backend is hidden, as CI installs without it, so that each refusal is
## shown to come before the binary check and not only on machines that lack it.
## Returns the function that puts it back.
hide_backend <- function() {
  no_bin <- tempfile("no_progsf90_")
  dir.create(no_bin)
  old_bin <- breedR.getOption("breedR.bin")
  breedR.setOption("breedR.bin", no_bin)
  function() {
    breedR.setOption("breedR.bin", old_bin)
    unlink(no_bin, recursive = TRUE)
  }
}

## A RENUMF90 output directory, as far as the REML entry points read it.
## par_content is what parse_renumf90() captured when renumf90() ran; the file
## is what the last fit from this object left behind, which is not the same
## thing (see the test on a rewritten parameter file below).
fake_renum <- function(par = character(0), file_par = par) {
  d <- tempfile("fake_renum_")
  dir.create(d)
  par_file <- file.path(d, "renf90.par")
  writeLines(c("DATAFILE", "renf90.dat", file_par), par_file)
  list(par_file = par_file, dir = d,
       par_content = c("DATAFILE", "renf90.dat", par))
}

class_opts <- hetres_options(group_col = 2, n_groups = 2, var_file = "hv")

test_that("remlf90() refuses class-based heterogeneous residuals", {
  restore <- hide_backend()
  on.exit(restore(), add = TRUE)
  expect_false(check_progsf90(quiet = TRUE))

  fit <- function(opts)
    remlf90(phe_X ~ gg, data = globulus, progsf90.options = opts)

  msg <- tryCatch(fit(class_opts), error = conditionMessage)
  expect_match(msg, "hetres_int", fixed = TRUE)
  expect_match(msg, "gibbsf90(hetres_int", fixed = TRUE)

  ## however the option is written
  expect_error(fit("hetres_int 2 2"), "hetres_int")
  expect_error(fit(c("sol se", "  hetres_int 2 2")), "hetres_int")

  ## the covariate-based form is estimated by AI-REML and is not refused:
  ## it goes on to the binary check
  expect_error(fit(hetres_options(covariate_cols = 3, initial = c(1, 0.01))),
               "Binary dependencies missing")
  ## and the match is on the whole option name
  expect_error(fit("hetres_integer 2 2"), "Binary dependencies missing")
})

test_that("remlf90_from_renum() refuses class-based heterogeneous residuals", {
  restore <- hide_backend()
  on.exit(restore(), add = TRUE)

  msg <- tryCatch(remlf90_from_renum(fake_renum(),
                                     progsf90.options = class_opts),
                  error = conditionMessage)
  expect_match(msg, "hetres_int", fixed = TRUE)

  ## RENUMF90 copies the OPTION lines it is given into renf90.par
  expect_error(remlf90_from_renum(fake_renum("OPTION hetres_int 2 2")),
               "hetres_int")

  ## control: no class-based option, so on to the binary check
  expect_error(remlf90_from_renum(fake_renum()), "not installed")

  ## A fit from a renum object rewrites its renf90.par in place, so a previous
  ## gibbsf90_from_renum() with class-based residuals leaves OPTION hetres_int
  ## in the file. That is the Gibbs fit's option, not this call's: refusing on
  ## it would tell a user who asked for a plain REML fit to go and use
  ## gibbsf90(), which is exactly what they just did.
  expect_error(remlf90_from_renum(fake_renum(file_par = "OPTION hetres_int 2 2")),
               "not installed")
})

test_that("validate_prediction() refuses class-based heterogeneous residuals", {
  restore <- hide_backend()
  on.exit(restore(), add = TRUE)

  expect_error(validate_prediction(fake_renum(), validation_ids = 1:2,
                                   progsf90.options = class_opts),
               "hetres_int")
  expect_error(validate_prediction(fake_renum("OPTION hetres_int 2 2"),
                                   validation_ids = 1:2),
               "hetres_int")

  ## as above: an option another fit left in the file is not this call's
  msg <- tryCatch(validate_prediction(fake_renum(file_par = "OPTION hetres_int 2 2"),
                                      validation_ids = 1:2),
                  error = conditionMessage)
  expect_false(grepl("hetres_int", msg, fixed = TRUE))
})
