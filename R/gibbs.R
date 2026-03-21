### GIBBSF90+ — Bayesian variance component estimation via Gibbs sampling ###


#' Bayesian variance component estimation via Gibbs sampling
#'
#' Fits a linear mixed model using Gibbs sampling via GIBBSF90+. Supports
#' threshold/categorical traits, heterogeneous residual variances, and
#' genomic models (ssGBLUP). Provides posterior distributions for variance
#' components and breeding values.
#'
#' @inheritParams remlf90
#' @param n_samples integer. Total number of MCMC samples (default 10000).
#' @param burnin integer. Number of initial samples to discard (default 0).
#'   Set to 0 for a first run, then inspect convergence with
#'   \code{\link{postgibbsf90}} to determine an appropriate burn-in.
#' @param thin integer. Save every n-th sample (default 10).
#' @param cat integer vector or NULL. Category specification per trait.
#'   0 = linear, 2 = binary, 3+ = ordinal with that many categories.
#'   Phenotypes for categorical traits must be numbered from 1 (0 = missing).
#' @param thresholds numeric vector or NULL. Fixed thresholds for ordinal
#'   traits. Not needed for binary traits.
#' @param prior numeric vector or NULL. Degree of belief for each effect.
#'   Format: eff1 db1 eff2 db2 ... -1 db_residual.
#' @param seed integer vector of length 2 or NULL. Random seeds.
#' @param fixed_var character or NULL. "mean" or "all" — estimate solutions
#'   with known (fixed) variances from the parameter file.
#' @param save_samples character or NULL. "mean" (posterior means + SDs only)
#'   or "all" (store all samples — creates large files).
#' @param save_effects integer vector or NULL. Which effects to save.
#' @param cont logical. Continue from a previous run (default FALSE).
#' @param save_halfway integer or NULL. Save checkpoint every N rounds.
#' @param hetres_int list or NULL. Heterogeneous residual variances.
#'   List with \code{col} and \code{n}.
#' @param residual_var numeric or NULL. Fixed residual variance for
#'   categorical traits.
#' @return A list with components:
#'   \describe{
#'     \item{solutions}{data.frame of posterior means + SDs for solutions}
#'     \item{samples}{matrix of variance component samples per round}
#'     \item{deviance}{numeric vector of -2logL per round}
#'     \item{output}{character vector of program stdout}
#'     \item{dir}{path to working directory with all output files}
#'     \item{n_samples}{total samples requested}
#'     \item{burnin}{burn-in used}
#'     \item{thin}{thinning interval used}
#'   }
#' @examples
#' \dontrun{
#' # Linear model with Gibbs sampling
#' res <- gibbsf90(phe_X ~ gg,
#'   genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),
#'   data = dat, n_samples = 50000, burnin = 5000, thin = 10)
#'
#' # Posterior means for variance components
#' colMeans(res$samples)
#'
#' # Threshold model for binary survival
#' res_bin <- gibbsf90(survival ~ site,
#'   genetic = list(model = 'add_animal', pedigree = ped, id = 'self'),
#'   data = dat, n_samples = 100000, burnin = 20000, thin = 10,
#'   cat = 2)
#' }
#' @export
gibbsf90 <- function(fixed,
                      random = NULL,
                      genetic = NULL,
                      spatial = NULL,
                      generic = NULL,
                      genomic = NULL,
                      data,
                      var.ini = NULL,
                      n_samples = 10000L,
                      burnin = 0L,
                      thin = 10L,
                      cat = NULL,
                      thresholds = NULL,
                      prior = NULL,
                      seed = NULL,
                      fixed_var = NULL,
                      save_samples = "mean",
                      save_effects = NULL,
                      cont = FALSE,
                      save_halfway = NULL,
                      hetres_int = NULL,
                      residual_var = NULL,
                      breedR.bin = breedR.getOption("breedR.bin"),
                      progsf90.options = NULL,
                      weights = NULL) {

  mc <- match.call()

  # Validate MCMC parameters
  n_samples <- as.integer(n_samples)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  if (n_samples < 1) stop("'n_samples' must be positive.", call. = FALSE)
  if (burnin < 0) stop("'burnin' must be non-negative.", call. = FALSE)
  if (thin < 1) stop("'thin' must be positive.", call. = FALSE)

  # Check binary
  gibbs_name <- "gibbsf90+"
  if (breedR.os.type() == 'windows') gibbs_name <- paste0(gibbs_name, ".exe")
  if (!file.exists(file.path(breedR.bin, gibbs_name)))
    stop("gibbsf90+ binary not found. Run install_genomic_programs().",
         call. = FALSE)

  ## --- Shared model-building pipeline (same as remlf90) ---

  method <- 'em'  # placeholder — Gibbs doesn't use REML method

  # Build model frame
  mf <- build.mf(mc)
  responsem <- as.matrix(stats::model.response(mf))

  # Validate components
  if (!is.null(genetic))
    genetic <- do.call('check_genetic',
      c(genetic, list(data = data, response = responsem)))

  if (!is.null(spatial))
    spatial <- do.call('check_spatial',
      c(spatial, list(data = data, response = responsem)))

  if (!is.null(generic))
    generic <- check_generic(generic, response = responsem)

  if (!is.null(genomic)) {
    if (is.null(genetic))
      stop("'genomic' requires a 'genetic' component.", call. = FALSE)
    genomic <- check_genomic(genomic)
  }

  var.ini <- check_var.ini(var.ini, random, responsem)

  # Build effects
  effects <- build.effects(mf, genetic, spatial, generic, var.ini)

  ## --- Gibbs-specific options ---

  gibbs_opts <- build_gibbs_options(
    cat = cat, thresholds = thresholds, prior = prior, seed = seed,
    fixed_var = fixed_var, save_samples = save_samples,
    save_effects = save_effects, cont = cont,
    save_halfway = save_halfway, hetres_int = hetres_int,
    residual_var = residual_var
  )

  genomic_opts <- if (!is.null(genomic)) build_genomic_options(genomic)
                  else NULL

  # Build progsf90 object (no 'sol se' or 'method VCE' — Gibbs doesn't use these)
  pf90 <- progsf90(mf, weights = weights, effects,
                    opt = c(gibbs_opts, genomic_opts, progsf90.options),
                    res.var.ini = var.ini$residuals)

  # Write files
  tmpdir <- tempdir()
  write.progsf90(pf90, dir = tmpdir)

  # Genomic preprocessing (if needed)
  pregs_out <- NULL
  if (!is.null(genomic)) {
    pipeline <- run_pregsf90_pipeline(genomic, genomic_opts, pf90,
                                      tmpdir, breedR.bin)
    pf90 <- pipeline$pf90
    pregs_out <- pipeline$pregs_out
    write.progsf90(pf90, dir = tmpdir)
  }

  ## --- Run GIBBSF90+ ---

  gibbs_src <- file.path(breedR.bin, gibbs_name)
  gibbs_bin <- file.path(tmpdir, gibbs_name)
  file.copy(gibbs_src, gibbs_bin, overwrite = TRUE)

  cdir <- setwd(tmpdir)
  on.exit({
    setwd(cdir)
    unlink(gibbs_bin)
  })

  args <- c("parameters",
            "--samples", as.character(n_samples),
            "--burnin", as.character(burnin),
            "--interval", as.character(thin))

  out <- system2(file.path(".", gibbs_name), args = args,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(tmpdir, "gibbsf90.out"))
    stop("gibbsf90+ failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(out, file.path(tmpdir, "gibbsf90.out"))

  ## --- Parse results ---

  result <- parse_gibbs_results(tmpdir)
  result$output <- out
  result$call <- mc
  result$n_samples <- n_samples
  result$burnin <- burnin
  result$thin <- thin
  result$dir <- tmpdir

  if (!is.null(genomic)) {
    result$genomic <- parse_pregsf90_qc(tmpdir)
  }

  return(result)
}


## Parse GIBBSF90+ output files
parse_gibbs_results <- function(dir) {

  result <- list()

  # Final solutions (posterior means + SDs)
  # Format: header line, then "trait effect level solution SD"
  sol_file <- file.path(dir, "final_solutions")
  if (file.exists(sol_file) && file.info(sol_file)$size > 0) {
    sol <- tryCatch(
      utils::read.table(sol_file, header = TRUE, row.names = NULL),
      error = function(e) NULL
    )
    if (!is.null(sol)) {
      col_names <- c("trait", "effect", "level", "solution", "sd")
      names(sol) <- col_names[seq_len(ncol(sol))]
      result$solutions <- sol
    }
  }

  # Gibbs samples (variance components per round)
  # Format: 3 header lines (effect/trait mapping), then alternating pairs:
  #   round_number  n_components
  #   value1  value2  ...  valueN
  samples_file <- file.path(dir, "gibbs_samples")
  if (file.exists(samples_file) && file.info(samples_file)$size > 0) {
    raw_lines <- readLines(samples_file)
    if (length(raw_lines) > 3) {
      data_lines <- raw_lines[-(1:3)]  # skip 3 header lines
      # Odd lines = "round_number  n_components"
      # Even lines = "value1  value2  ..."
      value_idx <- seq(2, length(data_lines), by = 2)
      round_idx <- seq(1, length(data_lines), by = 2)
      if (length(value_idx) > 0) {
        result$samples <- tryCatch({
          con <- textConnection(data_lines[value_idx])
          mat <- as.matrix(utils::read.table(con))
          close(con)
          round_nums <- as.integer(
            sub("\\s+.*", "", trimws(data_lines[round_idx])))
          rownames(mat) <- round_nums[seq_len(nrow(mat))]
          mat
        }, error = function(e) NULL)
      }
    }
  }

  # Deviance per round
  deviance_file <- file.path(dir, "fort.99")
  if (file.exists(deviance_file) && file.info(deviance_file)$size > 0) {
    result$deviance <- tryCatch(
      scan(deviance_file, quiet = TRUE),
      error = function(e) NULL
    )
  }

  return(result)
}


#' Build OPTION lines for GIBBSF90+
#'
#' Converts Gibbs-specific R arguments into OPTION strings for the
#' GIBBSF90+ parameter file.
#'
#' @param cat integer vector or NULL. Category specification per trait.
#'   0 = linear, 2+ = categorical with that many categories.
#' @param thresholds numeric vector or NULL. Fixed thresholds for
#'   categorical traits. Not needed for binary (2-category) traits.
#' @param prior numeric vector or NULL. Degree of belief for each random
#'   effect and residual. Format: eff1 db1 eff2 db2 ... -1 db_residual.
#' @param seed integer vector of length 2 or NULL. Seeds for reproducibility.
#' @param fixed_var character or NULL. If "mean", compute posterior means
#'   with known (fixed) variances. If "all", also store all samples.
#' @param save_samples character or NULL. "mean" = posterior means + SDs
#'   only. "all" = store all samples (large files).
#' @param save_effects integer vector or NULL. Which effects to save
#'   solutions for. NULL = all.
#' @param cont logical. Continue from a previous run (requires checkpoint
#'   files from save_halfway).
#' @param save_halfway integer or NULL. Save checkpoint every N rounds
#'   for cold restart.
#' @param hetres_int list or NULL. Heterogeneous residual variances by
#'   class. List with \code{col} (data column) and \code{n} (number of
#'   classes).
#' @param residual_var numeric or NULL. Fixed residual variance for
#'   categorical traits (automatically 1 for binary).
#' @return Character vector of OPTION strings (without "OPTION" prefix).
build_gibbs_options <- function(cat = NULL,
                                 thresholds = NULL,
                                 prior = NULL,
                                 seed = NULL,
                                 fixed_var = NULL,
                                 save_samples = NULL,
                                 save_effects = NULL,
                                 cont = FALSE,
                                 save_halfway = NULL,
                                 hetres_int = NULL,
                                 residual_var = NULL) {

  opts <- character(0)

  # Categorical trait specification
  if (!is.null(cat)) {
    if (!is.numeric(cat) || any(cat < 0))
      stop("'cat' must be a non-negative integer vector.", call. = FALSE)
    opts <- c(opts, paste("cat", paste(as.integer(cat), collapse = " ")))
  }

  # Fixed thresholds
  if (!is.null(thresholds)) {
    if (!is.numeric(thresholds))
      stop("'thresholds' must be numeric.", call. = FALSE)
    opts <- c(opts, paste("thresholds",
                           paste(format(thresholds, scientific = FALSE),
                                 collapse = " ")))
  }

  # Prior degrees of belief
  if (!is.null(prior)) {
    if (!is.numeric(prior))
      stop("'prior' must be a numeric vector.", call. = FALSE)
    opts <- c(opts, paste("prior",
                           paste(prior, collapse = " ")))
  }

  # Random seeds
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 2)
      stop("'seed' must be an integer vector of length 2.", call. = FALSE)
    opts <- c(opts, paste("seed", seed[1], seed[2]))
  }

  # Fixed variance mode (known variances, estimate solutions only)
  if (!is.null(fixed_var)) {
    if (!fixed_var %in% c("mean", "all"))
      stop("'fixed_var' must be 'mean' or 'all'.", call. = FALSE)
    fv_line <- paste("fixed_var", fixed_var)
    if (!is.null(save_effects))
      fv_line <- paste(fv_line, paste(as.integer(save_effects), collapse = " "))
    opts <- c(opts, fv_line)
  } else if (!is.null(save_samples)) {
    # Solution saving (when variances are estimated)
    if (!save_samples %in% c("mean", "all"))
      stop("'save_samples' must be 'mean' or 'all'.", call. = FALSE)
    sol_line <- paste("solution", save_samples)
    if (!is.null(save_effects))
      sol_line <- paste(sol_line, paste(as.integer(save_effects), collapse = " "))
    opts <- c(opts, sol_line)
  }

  # Continue from previous run
  if (isTRUE(cont))
    opts <- c(opts, "cont 1")

  # Save checkpoint for cold restart
  if (!is.null(save_halfway)) {
    if (!is.numeric(save_halfway) || save_halfway <= 0)
      stop("'save_halfway' must be a positive integer.", call. = FALSE)
    opts <- c(opts, paste("save_halfway_samples", as.integer(save_halfway)))
  }

  # Heterogeneous residual variances (native in Gibbs)
  if (!is.null(hetres_int)) {
    if (!is.list(hetres_int) || is.null(hetres_int$col) || is.null(hetres_int$n))
      stop("'hetres_int' must be a list with 'col' and 'n'.", call. = FALSE)
    opts <- c(opts, paste("hetres_int", hetres_int$col, hetres_int$n))
  }

  # Fixed residual variance for categorical traits
  if (!is.null(residual_var))
    opts <- c(opts, paste("residual", residual_var))

  return(opts)
}


### POSTGIBBSF90 — Post-processing of Gibbs samples ###


#' Post-process Gibbs samples for convergence diagnostics
#'
#' Runs POSTGIBBSF90 on the output of \code{\link{gibbsf90}} to compute
#' posterior means, standard deviations, HPD intervals, effective sample
#' sizes, Geweke convergence diagnostics, and autocorrelations.
#'
#' Must be called in the same R session as \code{\link{gibbsf90}} since
#' it reads files from \code{tempdir()}.
#'
#' @param gibbs_result output from \code{\link{gibbsf90}}, or a character
#'   path to the directory containing \code{gibbs_samples} and
#'   \code{fort.99}.
#' @param burnin integer. Burn-in for post-processing. If NULL, defaults
#'   to the burn-in used in the \code{gibbsf90()} call.
#' @param thin integer. Thinning interval for post-processing. If NULL,
#'   defaults to the thinning used in the \code{gibbsf90()} call.
#' @return A list with:
#'   \describe{
#'     \item{posterior_means}{data.frame with posterior means per parameter}
#'     \item{posterior_sd}{data.frame with posterior SDs per parameter}
#'     \item{hpd}{data.frame with 95\% HPD intervals}
#'     \item{effective_size}{numeric vector of effective sample sizes}
#'     \item{geweke}{numeric vector of Geweke convergence diagnostics}
#'     \item{autocorrelations}{data.frame of autocorrelations at lag 1, 10, 50}
#'     \item{log_marginal_density}{numeric, log marginal density for Bayes factor}
#'     \item{diagnostics_raw}{character vector of raw postout file content}
#'     \item{output}{character vector of POSTGIBBSF90 stdout}
#'   }
#' @export
postgibbsf90 <- function(gibbs_result,
                          burnin = NULL,
                          thin = NULL) {

  # Determine working directory
  if (is.character(gibbs_result)) {
    dir <- gibbs_result
  } else if (is.list(gibbs_result) && !is.null(gibbs_result$dir)) {
    dir <- gibbs_result$dir
    if (is.null(burnin)) burnin <- gibbs_result$burnin
    if (is.null(thin)) thin <- gibbs_result$thin
  } else {
    stop("'gibbs_result' must be a gibbsf90() result or a directory path.",
         call. = FALSE)
  }

  if (is.null(burnin)) burnin <- 0L
  if (is.null(thin)) thin <- 1L

  # Check required input files
  if (!file.exists(file.path(dir, "gibbs_samples")))
    stop("gibbs_samples file not found in ", dir, call. = FALSE)
  if (!file.exists(file.path(dir, "fort.99")))
    stop("fort.99 file not found in ", dir, call. = FALSE)

  # Run POSTGIBBSF90
  out <- run_postgibbsf90_bin(dir, breedR.getOption('breedR.bin'),
                               burnin, thin)
  writeLines(out, file.path(dir, "postgibbsf90.out"))

  # Parse output
  result <- parse_postgibbsf90(dir, out)
  result$output <- out

  return(result)
}


## Execute POSTGIBBSF90 binary
## Interactive program — pipe parameter file name, burnin, thin, "0" via stdin
run_postgibbsf90_bin <- function(dir, bin_path, burnin, thin) {

  postg_name <- "postgibbsf90"
  if (breedR.os.type() == 'windows') postg_name <- paste0(postg_name, ".exe")
  postg_src <- file.path(bin_path, postg_name)

  if (!file.exists(postg_src))
    stop("postgibbsf90 binary not found. Run install_genomic_programs().",
         call. = FALSE)

  # Copy binary to working directory (DLL isolation)
  postg_bin <- file.path(dir, postg_name)
  file.copy(postg_src, postg_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(postg_bin)
  })

  # POSTGIBBSF90 prompts:
  # 1. "name of parameter file?" → parameter file name
  # 2. "Burn-in?" → integer
  # 3. "Give n to read every n-th sample?" → integer
  # 4. "Choose a graph..." → "0" to exit
  input_lines <- c("parameters",
                    as.character(as.integer(burnin)),
                    as.character(as.integer(thin)),
                    "0")

  out <- system2(file.path(".", postg_name), input = input_lines,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "postgibbsf90.out"))
    stop("postgibbsf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  return(out)
}


## Parse POSTGIBBSF90 output files and stdout
parse_postgibbsf90 <- function(dir, out) {

  result <- list()

  # Parse postmean file
  postmean_file <- file.path(dir, "postmean")
  if (file.exists(postmean_file) && file.info(postmean_file)$size > 0) {
    result$posterior_means <- tryCatch(
      utils::read.table(postmean_file),
      error = function(e) NULL)
  }

  # Parse postsd file
  postsd_file <- file.path(dir, "postsd")
  if (file.exists(postsd_file) && file.info(postsd_file)$size > 0) {
    result$posterior_sd <- tryCatch(
      utils::read.table(postsd_file),
      error = function(e) NULL)
  }

  # Parse postout file — the main diagnostics table
  postout_file <- file.path(dir, "postout")
  if (file.exists(postout_file) && file.info(postout_file)$size > 0) {
    postout_lines <- readLines(postout_file)
    result$diagnostics_raw <- postout_lines
    result <- c(result, parse_postout_tables(postout_lines))
  }

  # Parse postgibbs_samples (thinned samples for external use)
  pg_samples_file <- file.path(dir, "postgibbs_samples")
  if (file.exists(pg_samples_file) && file.info(pg_samples_file)$size > 0) {
    result$samples <- tryCatch(
      as.matrix(utils::read.table(pg_samples_file)),
      error = function(e) NULL)
  }

  # Extract log marginal density from stdout
  lmd_line <- grep("log\\(p\\)", out, value = TRUE)
  if (length(lmd_line) > 0) {
    lmd_val <- as.numeric(sub(".*=\\s*", "", trimws(lmd_line[1])))
    if (!is.na(lmd_val)) result$log_marginal_density <- lmd_val
  }

  return(result)
}


## Parse the postout diagnostics tables
## Two sections: "Monte Carlo Error" and "Posterior Standard Deviation"
parse_postout_tables <- function(lines) {

  result <- list()

  # Find numeric data lines (lines starting with position numbers)
  # Format: Pos eff1 eff2 trt1 trt2 value1 value2 ...
  data_lines <- grep("^\\s*[0-9]+\\s+[0-9]+", lines, value = TRUE)

  if (length(data_lines) == 0) return(result)

  # Split into two tables: MCE table and PSD table
  # They're separated by the "Posterior Standard Deviation" header
  psd_header <- grep("Posterior Standard Deviation", lines)

  if (length(psd_header) > 0) {
    # Find which data lines come before and after the PSD header
    data_line_nums <- grep("^\\s*[0-9]+\\s+[0-9]+", lines)
    mce_lines <- data_line_nums[data_line_nums < psd_header[1]]
    psd_lines <- data_line_nums[data_line_nums > psd_header[1]]

    # Parse MCE table
    if (length(mce_lines) > 0) {
      mce_data <- tryCatch({
        con <- textConnection(lines[mce_lines])
        df <- utils::read.table(con)
        close(con)
        # Columns: Pos, eff1, eff2, trt1, trt2, MCE, Mean, HPD_low, HPD_high,
        #          EffSampleSize, Median, Mode, IndChainSize
        col_names <- c("pos", "eff1", "eff2", "trt1", "trt2",
                        "mce", "mean", "hpd_lower", "hpd_upper",
                        "effective_size", "median", "mode",
                        "independent_chain_size")
        if (ncol(df) <= length(col_names))
          names(df) <- col_names[seq_len(ncol(df))]
        df
      }, error = function(e) NULL)

      if (!is.null(mce_data)) {
        result$mce_table <- mce_data
        result$hpd <- mce_data[, c("pos", "eff1", "eff2", "trt1", "trt2",
                                     "hpd_lower", "hpd_upper")]
        result$effective_size <- mce_data$effective_size
      }
    }

    # Parse PSD table
    if (length(psd_lines) > 0) {
      psd_data <- tryCatch({
        con <- textConnection(lines[psd_lines])
        df <- utils::read.table(con)
        close(con)
        # Columns: Pos, eff1, eff2, trt1, trt2, PSD, Mean, PSD_low, PSD_high,
        #          Geweke, AutoCorr_lag1, lag10, lag50, IndBatches
        col_names <- c("pos", "eff1", "eff2", "trt1", "trt2",
                        "psd", "mean", "psd_lower", "psd_upper",
                        "geweke", "autocorr_lag1", "autocorr_lag10",
                        "autocorr_lag50", "independent_batches")
        if (ncol(df) <= length(col_names))
          names(df) <- col_names[seq_len(ncol(df))]
        df
      }, error = function(e) NULL)

      if (!is.null(psd_data)) {
        result$psd_table <- psd_data
        result$geweke <- psd_data$geweke
        result$autocorrelations <- psd_data[, c("pos", "autocorr_lag1",
                                                  "autocorr_lag10",
                                                  "autocorr_lag50")]
      }
    }
  }

  return(result)
}


#' Fit a Gibbs sampling model using RENUMF90 output
#'
#' Takes the output of \code{\link{renumf90}} and runs GIBBSF90+ for
#' Bayesian variance component estimation. Supports threshold/categorical
#' traits via the \code{cat} parameter.
#'
#' @param renum output from \code{\link{renumf90}}.
#' @param n_samples integer. Total MCMC samples.
#' @param burnin integer. Burn-in samples to discard.
#' @param thin integer. Thinning interval.
#' @param cat integer vector or NULL. Category specification per trait
#'   (0 = linear, 2 = binary, 3+ = ordinal).
#' @param progsf90.options character vector. Additional OPTIONS.
#' @param prior numeric vector or NULL. Degree of belief.
#' @param seed integer vector of length 2 or NULL.
#' @param save_samples character. "mean" or "all".
#' @return A list with solutions, variance samples, deviance, and
#'   the RENUMF90 output.
#' @export
gibbsf90_from_renum <- function(renum,
                                 n_samples = 10000L,
                                 burnin = 0L,
                                 thin = 10L,
                                 cat = NULL,
                                 progsf90.options = NULL,
                                 prior = NULL,
                                 seed = NULL,
                                 save_samples = "mean") {

  if (is.null(renum$par_file) || !file.exists(renum$par_file))
    stop("RENUMF90 parameter file not found. Run renumf90() first.",
         call. = FALSE)

  dir <- renum$dir
  bin_path <- breedR.getOption('breedR.bin')

  gibbs_name <- "gibbsf90+"
  if (breedR.os.type() == 'windows') gibbs_name <- paste0(gibbs_name, ".exe")
  if (!file.exists(file.path(bin_path, gibbs_name)))
    stop("gibbsf90+ binary not found. Run install_genomic_programs().",
         call. = FALSE)

  # Build Gibbs options
  gibbs_opts <- build_gibbs_options(
    cat = cat, prior = prior, seed = seed, save_samples = save_samples
  )

  # Read RENUMF90 parameter file and add options
  par_lines <- readLines(renum$par_file)
  par_lines <- c(par_lines,
                  if (length(gibbs_opts) > 0) paste("OPTION", gibbs_opts),
                  if (!is.null(progsf90.options))
                    paste("OPTION", progsf90.options))
  writeLines(par_lines, file.path(dir, "renf90.par"))

  # Run GIBBSF90+
  gibbs_src <- file.path(bin_path, gibbs_name)
  gibbs_bin <- file.path(dir, gibbs_name)
  file.copy(gibbs_src, gibbs_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(gibbs_bin)
  })

  args <- c("renf90.par",
            "--samples", as.character(as.integer(n_samples)),
            "--burnin", as.character(as.integer(burnin)),
            "--interval", as.character(as.integer(thin)))

  out <- system2(file.path(".", gibbs_name), args = args,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status')))
    stop("gibbsf90+ failed with exit code ", attr(out, 'status'),
         ".\nOutput:\n", paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)

  writeLines(out, file.path(dir, "gibbsf90.out"))

  # Parse results
  result <- parse_gibbs_results(dir)
  result$output <- out
  result$n_samples <- as.integer(n_samples)
  result$burnin <- as.integer(burnin)
  result$thin <- as.integer(thin)
  result$dir <- dir
  result$renum <- renum

  return(result)
}
