### VALIDATIONF90 — Prediction validation ###


#' Validate genomic/genetic predictions
#'
#' Implements the LR validation method (Legarra & Reverter, 2018) by
#' comparing solutions from a full dataset against a partial (training)
#' dataset. Computes bias, dispersion, accuracy, and optionally
#' predictive ability and confidence intervals.
#'
#' The validation workflow is:
#' \enumerate{
#'   \item Fit the model with the full dataset (all animals)
#'   \item Fit the model with a partial dataset (remove validation animals'
#'     phenotypes)
#'   \item Call \code{validationf90()} with both solution sets and a list
#'     of validation animal IDs
#' }
#'
#' @param par_file character. Path to the BLUPF90+ parameter file
#'   (e.g., \code{renf90.par} or the \code{parameters} file from
#'   \code{remlf90()}).
#' @param solutions_whole character. Path to the solutions file from the
#'   full dataset analysis.
#' @param solutions_partial character. Path to the solutions file from the
#'   partial (training) dataset analysis.
#' @param validation_ids integer or character vector. Renumbered IDs of
#'   validation animals, or path to a file containing them (one per line).
#' @param effect integer. The effect number to validate (typically the
#'   genetic effect).
#' @param predictive_ability logical. If TRUE, also compute predictive
#'   ability (correlation between predictions and adjusted phenotypes).
#'   Requires \code{yhat_residual} file from \code{predictf90}
#'   (default FALSE).
#' @param se character or NULL. Standard error method: \code{"exact"}
#'   (needs PEV matrices), \code{"approx"} (needs reliabilities),
#'   \code{"boot"} (bootstrap), or NULL (no SEs, default).
#' @param focal_var numeric vector or NULL. Genetic variance in the
#'   validation set for each trait. If NULL (default), estimated from
#'   inbreeding.
#' @param dir character. Working directory containing input files.
#'   Default \code{tempdir()}.
#' @return A list with validation statistics:
#'   \describe{
#'     \item{output}{character vector of raw program output}
#'     \item{statistics}{parsed validation statistics (if available)}
#'   }
#' @export
validationf90 <- function(par_file,
                           solutions_whole,
                           solutions_partial,
                           validation_ids,
                           effect,
                           predictive_ability = FALSE,
                           se = NULL,
                           focal_var = NULL,
                           dir = tempdir()) {

  # Validate inputs
  for (f in c(par_file, solutions_whole, solutions_partial)) {
    if (!file.exists(f))
      stop("File not found: ", f, call. = FALSE)
  }

  bin_path <- breedR.getOption('breedR.bin')
  val_name <- "validationf90"
  if (breedR.os.type() == 'windows') val_name <- paste0(val_name, ".exe")
  val_src <- file.path(bin_path, val_name)

  if (!file.exists(val_src))
    stop("validationf90 binary not found. Run install_genomic_programs().",
         call. = FALSE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Copy parameter file
  par_dest <- file.path(dir, "renf90.par")
  file.copy(par_file, par_dest, overwrite = TRUE)

  # Copy solutions files with the expected names
  file.copy(solutions_whole,
            file.path(dir, "solutions_whole"), overwrite = TRUE)
  file.copy(solutions_partial,
            file.path(dir, "solutions_partial"), overwrite = TRUE)

  # Write validation ID file
  if (is.character(validation_ids) && length(validation_ids) == 1 &&
      file.exists(validation_ids)) {
    # It's a file path
    file.copy(validation_ids,
              file.path(dir, "validation_ids.txt"), overwrite = TRUE)
  } else {
    # It's a vector of IDs
    writeLines(as.character(validation_ids),
               file.path(dir, "validation_ids.txt"))
  }

  # Add validation OPTIONs to parameter file
  par_lines <- readLines(par_dest)
  # Remove any existing validation options
  par_lines <- par_lines[!grepl("^OPTION validation|^OPTION predictive|^OPTION focal_var|^OPTION se ",
                                 par_lines)]

  # Add validation option
  par_lines <- c(par_lines,
                  paste("OPTION validation", effect, "validation_ids.txt"))

  if (isTRUE(predictive_ability))
    par_lines <- c(par_lines, "OPTION predictive_ability")

  if (!is.null(se))
    par_lines <- c(par_lines, paste("OPTION se", se))

  if (!is.null(focal_var))
    par_lines <- c(par_lines,
                    paste("OPTION focal_var usr",
                          paste(focal_var, collapse = " ")))

  writeLines(par_lines, par_dest)

  # Copy binary to working directory (DLL isolation)
  val_bin <- file.path(dir, val_name)
  file.copy(val_src, val_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(val_bin)
  })

  out <- system2(file.path(".", val_name), input = "renf90.par",
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "validationf90.out"))
    stop("validationf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(out, file.path(dir, "validationf90.out"))

  # Parse output
  result <- list(output = out)

  # Extract key statistics from output text
  result$statistics <- parse_validation_output(out)

  return(result)
}


#' Run a full validation analysis
#'
#' Convenience function that fits the model twice (full and partial data),
#' then runs validationf90. This automates the multi-step validation workflow.
#'
#' @param renum output from \code{\link{renumf90}}.
#' @param validation_ids integer vector. Renumbered IDs of validation animals.
#' @param effect integer. Effect number to validate.
#' @param method either 'ai' or 'em'.
#' @param progsf90.options character vector. Additional options.
#' @param partial_data_file character. Path to the partial dataset (with
#'   validation animals' phenotypes removed). If NULL, the function creates
#'   it by setting validation animals' observations to missing (0).
#' @return Validation statistics from \code{\link{validationf90}}.
#' @export
validate_prediction <- function(renum,
                                 validation_ids,
                                 effect = 2L,
                                 method = c('ai', 'em'),
                                 progsf90.options = NULL,
                                 partial_data_file = NULL) {

  method <- match.arg(method)

  if (is.null(renum$par_file) || !file.exists(renum$par_file))
    stop("RENUMF90 output required. Run renumf90() first.", call. = FALSE)

  dir <- renum$dir
  bin_path <- breedR.getOption('breedR.bin')
  blupf90_name <- progsf90_files(breedR.os.type())

  # Build method options
  method_opts <- 'method VCE'
  if (method == 'em') method_opts <- c(method_opts, 'EM-REML')
  extra_opts <- c("sol se", method_opts, progsf90.options)

  # Read base parameter file
  par_lines <- readLines(renum$par_file)
  par_with_opts <- c(par_lines, paste("OPTION", extra_opts))

  # --- Step 1: Fit with WHOLE data ---
  message("Fitting model with full data...")
  whole_dir <- file.path(dir, "validation_whole")
  dir.create(whole_dir, showWarnings = FALSE, recursive = TRUE)

  # Copy all RENUMF90 files
  renum_files <- list.files(dir, full.names = TRUE)
  renum_files <- renum_files[!file.info(renum_files)$isdir]
  file.copy(renum_files, whole_dir, overwrite = TRUE)

  writeLines(par_with_opts, file.path(whole_dir, "renf90.par"))

  run_blupf90_in_dir(whole_dir, bin_path, blupf90_name)

  sol_whole <- file.path(whole_dir, "solutions")
  if (!file.exists(sol_whole))
    stop("BLUPF90+ did not produce solutions for full data.", call. = FALSE)

  # --- Step 2: Create partial data and fit ---
  message("Fitting model with partial (training) data...")
  partial_dir <- file.path(dir, "validation_partial")
  dir.create(partial_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(renum_files, partial_dir, overwrite = TRUE)

  if (is.null(partial_data_file)) {
    # Create partial data by setting validation animals' phenotypes to missing
    data_file <- file.path(partial_dir, "renf90.dat")
    dat <- readLines(data_file)

    # Parse the parameter file to find trait and effect column positions
    trait_line <- par_lines[grep("^OBSERVATION", par_lines) + 1]
    trait_cols <- as.integer(strsplit(trimws(trait_line), "\\s+")[[1]])

    # Find the data column for the genetic effect (the effect being validated).
    # The EFFECTS block lists one effect per line after the header.
    # Each line starts with the data column position(s).
    effects_start <- grep("^EFFECTS:", par_lines)
    if (length(effects_start) == 0)
      stop("Could not find EFFECTS block in parameter file.", call. = FALSE)

    # Read the n_effects lines after the EFFECTS header
    n_effects_line <- par_lines[grep("^NUMBER_OF_EFFECTS", par_lines) + 1]
    n_effects <- as.integer(trimws(n_effects_line))

    effect_lines <- par_lines[(effects_start + 1):(effects_start + n_effects)]
    # The validated effect's line gives us the data column position
    # (first number on the line for single-trait, or the first n_traits numbers)
    effect_fields <- strsplit(trimws(effect_lines[effect]), "\\s+")[[1]]
    # First field is the data column position for trait 1
    genetic_col <- as.integer(effect_fields[1])

    if (is.na(genetic_col) || genetic_col == 0)
      stop("Could not determine data column for effect ", effect,
           ". Check parameter file.", call. = FALSE)

    # Only check the genetic effect column for validation IDs
    val_ids_int <- as.integer(validation_ids)
    dat_split <- strsplit(dat, "\\s+")
    n_zeroed <- 0L
    for (i in seq_along(dat_split)) {
      fields <- dat_split[[i]]
      fields <- fields[nchar(fields) > 0]
      animal_id <- as.integer(fields[genetic_col])
      if (!is.na(animal_id) && animal_id %in% val_ids_int) {
        for (tc in trait_cols) {
          fields[tc] <- "0"
        }
        n_zeroed <- n_zeroed + 1L
      }
      dat_split[[i]] <- fields
    }
    writeLines(vapply(dat_split, paste, "", collapse = " "), data_file)
    message("  Zeroed phenotypes for ", n_zeroed, " validation animals")
  } else {
    file.copy(partial_data_file,
              file.path(partial_dir, "renf90.dat"), overwrite = TRUE)
  }

  writeLines(par_with_opts, file.path(partial_dir, "renf90.par"))

  run_blupf90_in_dir(partial_dir, bin_path, blupf90_name)

  sol_partial <- file.path(partial_dir, "solutions")
  if (!file.exists(sol_partial))
    stop("BLUPF90+ did not produce solutions for partial data.",
         call. = FALSE)

  # --- Step 3: Run validationf90 ---
  message("Running validation...")
  validationf90(
    par_file = file.path(whole_dir, "renf90.par"),
    solutions_whole = sol_whole,
    solutions_partial = sol_partial,
    validation_ids = validation_ids,
    effect = effect,
    dir = dir
  )
}


## Run BLUPF90+ in a specific directory
run_blupf90_in_dir <- function(dir, bin_path, blupf90_name) {
  blupf90_src <- file.path(bin_path, blupf90_name)
  blupf90_bin <- file.path(dir, blupf90_name)
  file.copy(blupf90_src, blupf90_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(blupf90_bin)
  })

  out <- system2(file.path(".", blupf90_name), input = "renf90.par",
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status')))
    stop("BLUPF90+ failed in ", dir, ".\nOutput:\n",
         paste(utils::tail(out, 15), collapse = "\n"), call. = FALSE)

  writeLines(out, file.path(dir, "blupf90.out"))
  invisible(out)
}


## Parse validationf90 output text for key statistics
parse_validation_output <- function(out) {
  stats <- list()

  # Look for common validation statistics in the output
  # Bias (b0): intercept of regression of whole on partial
  b0_line <- grep("bias|b0|intercept", out, value = TRUE, ignore.case = TRUE)
  if (length(b0_line) > 0) stats$bias_lines <- b0_line

  # Dispersion (b1): slope
  b1_line <- grep("dispersion|slope|b1", out, value = TRUE, ignore.case = TRUE)
  if (length(b1_line) > 0) stats$dispersion_lines <- b1_line

  # Accuracy / correlation
  acc_line <- grep("accuracy|correlation|rho", out, value = TRUE,
                    ignore.case = TRUE)
  if (length(acc_line) > 0) stats$accuracy_lines <- acc_line

  # Predictive ability
  pa_line <- grep("predictive.ability", out, value = TRUE, ignore.case = TRUE)
  if (length(pa_line) > 0) stats$predictive_ability_lines <- pa_line

  # Store full output for manual inspection
  stats$raw <- out

  return(stats)
}
