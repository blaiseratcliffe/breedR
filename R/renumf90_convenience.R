### Convenience wrapper for RENUMF90 using R data.frames ###


#' Prepare data for BLUPF90 using R objects
#'
#' Converts R data.frames into the file format RENUMF90 expects,
#' maps column names to positions, detects factor vs numeric types,
#' and runs RENUMF90. This is a convenience layer over
#' \code{\link{renumf90}} that avoids manual column position
#' specification.
#'
#' For models with different effects per trait, nested covariates,
#' maternal effects, or random regressions, use \code{\link{renumf90}}
#' directly.
#'
#' @param data data.frame. The phenotypic data.
#' @param traits character vector. Column names of trait(s) in \code{data}.
#'   Must be numeric columns.
#' @param fixed character vector or NULL. Column names of fixed effects.
#'   Character and factor columns are treated as cross-classified.
#'   Numeric columns are treated as covariates unless listed in
#'   \code{factors}.
#' @param random character vector or NULL. Column names of random
#'   (diagonal) effects (e.g., block). Must be discrete grouping
#'   variables (factor, character, or integer codes).
#' @param pedigree data.frame or NULL. Three columns: animal, sire, dam.
#'   The first column name must exist in \code{data} — it identifies
#'   the genetic effect. Use 0 or NA for unknown parents.
#' @param factors character vector or NULL. Column names that should be
#'   treated as cross-classified factors even if numeric. Useful for
#'   numeric block or site codes.
#' @param covariates character vector or NULL. Column names that should
#'   be treated as covariates even if they could be factors. Overrides
#'   auto-detection.
#' @param genetic_variance numeric. Starting genetic (co)variance.
#'   Scalar (applied to diagonal), vector (diagonal elements), or
#'   full matrix. If NULL, defaults to 1/3 of phenotypic variance.
#' @param residual_variance numeric. Starting residual (co)variance.
#'   Same format as \code{genetic_variance}.
#'   If NULL, defaults to 2/3 of phenotypic variance.
#' @param random_variance named list or NULL. Starting variances for
#'   diagonal random effects. Names must match \code{random} column
#'   names. If NULL, defaults to 1/10 of phenotypic variance.
#' @param weights character or NULL. Column name for observation weights.
#' @param snp_file character or NULL. Path to genotype file.
#' @param inbreeding character or NULL. Inbreeding method (e.g.,
#'   "pedigree", "no-inbreeding").
#' @param ped_depth integer or NULL. Pedigree search depth (0 = all).
#' @param options character vector or NULL. Additional OPTION lines for
#'   RENUMF90.
#' @param dir character. Working directory.
#' @return Output from \code{\link{renumf90}} — a list with par_file,
#'   data_file, pedigree, inbreeding, etc.
#'
#' @examples
#' \dontrun{
#' # Single trait
#' renum <- renumf90_from_data(
#'   data = dat,
#'   traits = "height",
#'   fixed = c("site", "age"),
#'   pedigree = dat[, c("tree_id", "sire", "dam")],
#'   genetic_variance = 5,
#'   residual_variance = 10
#' )
#' res <- remlf90_from_renum(renum)
#'
#' # Multi-trait
#' renum <- renumf90_from_data(
#'   data = dat,
#'   traits = c("height", "diameter"),
#'   fixed = c("site"),
#'   random = c("block"),
#'   pedigree = dat[, c("id", "sire", "dam")],
#'   factors = c("block"),
#'   genetic_variance = matrix(c(5, 2, 2, 3), 2, 2),
#'   residual_variance = matrix(c(10, 3, 3, 8), 2, 2)
#' )
#' }
#' @export
renumf90_from_data <- function(data,
                                traits,
                                fixed = NULL,
                                random = NULL,
                                pedigree = NULL,
                                factors = NULL,
                                covariates = NULL,
                                genetic_variance = NULL,
                                residual_variance = NULL,
                                random_variance = NULL,
                                weights = NULL,
                                snp_file = NULL,
                                inbreeding = NULL,
                                ped_depth = NULL,
                                options = NULL,
                                dir = tempdir()) {

  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)

  n_traits <- length(traits)
  if (n_traits == 0)
    stop("'traits' must have at least one element.", call. = FALSE)

  # Validate column names exist in data
  all_cols <- c(traits, fixed, random, weights)
  missing_cols <- all_cols[!all_cols %in% names(data)]
  if (length(missing_cols) > 0)
    stop("Columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  # Validate trait columns are numeric
  for (t in traits) {
    if (!is.numeric(data[[t]]))
      stop("Trait column '", t, "' must be numeric.", call. = FALSE)
  }

  # Validate random effects are discrete (not continuous)
  if (!is.null(random)) {
    for (rnd in random) {
      col <- data[[rnd]]
      if (is.numeric(col) && !rnd %in% factors) {
        n_unique <- length(unique(col[!is.na(col)]))
        if (n_unique > 100)
          warning("Random effect '", rnd, "' has ", n_unique,
                  " unique numeric values. Did you mean to list it in ",
                  "'factors' or use it as a fixed covariate?", call. = FALSE)
      }
    }
  }

  # Create working directory early (before writing any files)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Validate pedigree
  genetic_id <- NULL
  if (!is.null(pedigree)) {
    if (!is.data.frame(pedigree) || ncol(pedigree) < 3)
      stop("'pedigree' must be a data.frame with at least 3 columns.",
           call. = FALSE)
    genetic_id <- names(pedigree)[1]
    if (!genetic_id %in% names(data))
      stop("Pedigree animal column '", genetic_id,
           "' not found in data.", call. = FALSE)
  }

  # Compute trait variances once (used for auto defaults)
  trait_vars <- vapply(traits,
    function(t) stats::var(data[[t]], na.rm = TRUE), numeric(1))

  # Determine column types for each effect
  classify_col <- function(col_name) {
    if (!is.null(covariates) && col_name %in% covariates) return("cov")
    if (!is.null(factors) && col_name %in% factors) return("cross")
    col <- data[[col_name]]
    if (is.factor(col) || is.character(col)) return("cross")
    if (is.numeric(col)) return("cov")
    return("cross")  # default
  }

  # Build the output data file columns in a fixed order:
  # [genetic_id] [fixed...] [random...] [weights] [traits...]
  # Carefully avoid duplicates while preserving position integrity.
  out_cols <- character(0)
  if (!is.null(genetic_id)) out_cols <- c(out_cols, genetic_id)
  if (!is.null(fixed)) {
    new_fixed <- setdiff(fixed, out_cols)
    out_cols <- c(out_cols, new_fixed)
  }
  if (!is.null(random)) {
    new_random <- setdiff(random, out_cols)
    out_cols <- c(out_cols, new_random)
  }
  if (!is.null(weights)) {
    if (!weights %in% out_cols) out_cols <- c(out_cols, weights)
  }
  new_traits <- setdiff(traits, out_cols)
  out_cols <- c(out_cols, new_traits)

  # Verify no column was lost — all requested columns must be in out_cols
  for (col in c(traits, fixed, random, weights, genetic_id)) {
    if (!is.null(col) && !col %in% out_cols)
      stop("Internal error: column '", col, "' missing from output layout.",
           call. = FALSE)
  }

  # Build position map: column name -> position in output file
  pos_map <- stats::setNames(seq_along(out_cols), out_cols)

  # Trait positions
  trait_positions <- unname(pos_map[traits])

  # Weight position
  weight_pos <- if (!is.null(weights)) unname(pos_map[weights]) else NULL

  # Build effects list for renumf90()
  effects <- list()

  # Fixed effects
  if (!is.null(fixed)) {
    for (fx in fixed) {
      etype <- classify_col(fx)
      # Always use alpha for cross-classified — RENUMF90 handles both
      # numeric and alphanumeric values with alpha form.
      form <- if (etype == "cross") "alpha" else NULL

      eff <- list(
        pos = rep(unname(pos_map[fx]), n_traits),
        type = etype
      )
      if (etype == "cross") eff$form <- form
      effects <- c(effects, list(eff))
    }
  }

  # Genetic effect (animal model)
  if (!is.null(pedigree)) {
    # Write pedigree file — only first 3 columns (animal, sire, dam)
    ped_file <- file.path(dir, "pedigree_data.txt")
    ped_out <- pedigree[, 1:3, drop = FALSE]
    ped_out[is.na(ped_out)] <- 0
    # Replace spaces in IDs with underscores to prevent field splitting
    for (j in seq_len(ncol(ped_out))) {
      if (is.character(ped_out[[j]]) || is.factor(ped_out[[j]]))
        ped_out[[j]] <- gsub("\\s+", "_", as.character(ped_out[[j]]))
    }
    utils::write.table(ped_out, ped_file, row.names = FALSE,
                        col.names = FALSE, quote = FALSE)

    # Auto-compute genetic variance if not provided
    genetic_variance <- expand_variance(genetic_variance, n_traits,
                                         trait_vars / 3)

    gen_eff <- list(
      pos = rep(unname(pos_map[genetic_id]), n_traits),
      type = "cross",
      form = "alpha",
      random = "animal",
      file = basename(ped_file),
      covariances = genetic_variance
    )
    effects <- c(effects, list(gen_eff))
  }

  # Diagonal random effects
  if (!is.null(random)) {
    for (rnd in random) {
      # Auto-compute random variance
      if (!is.null(random_variance) && rnd %in% names(random_variance)) {
        rnd_var <- expand_variance(random_variance[[rnd]], n_traits,
                                    trait_vars / 10)
      } else {
        rnd_var <- expand_variance(NULL, n_traits, trait_vars / 10)
      }

      rnd_eff <- list(
        pos = rep(unname(pos_map[rnd]), n_traits),
        type = "cross",
        form = "alpha",
        random = "diagonal",
        covariances = rnd_var
      )
      effects <- c(effects, list(rnd_eff))
    }
  }

  # Residual variance
  residual_variance <- expand_variance(residual_variance, n_traits,
                                        trait_vars * 2 / 3)

  # Build output data
  out_data <- data[, out_cols, drop = FALSE]

  # Replace spaces in character/factor columns to prevent field splitting
  for (j in seq_along(out_cols)) {
    col <- out_data[[j]]
    if (is.character(col) || is.factor(col))
      out_data[[j]] <- gsub("\\s+", "_", as.character(col))
  }

  # Replace NA in trait columns with 0 (RENUMF90 missing code)
  for (t in traits) {
    out_data[[t]][is.na(out_data[[t]])] <- 0
  }

  # Check for NAs in effect columns
  effect_cols <- c(fixed, random, genetic_id)
  if (length(effect_cols) > 0) {
    na_rows <- apply(out_data[, effect_cols, drop = FALSE], 1,
                      function(r) any(is.na(r)))
    if (any(na_rows)) {
      warning(sum(na_rows), " rows have NA in effect columns and will be ",
              "dropped.", call. = FALSE)
      out_data <- out_data[!na_rows, , drop = FALSE]
    }
  }

  if (nrow(out_data) == 0)
    stop("No data remaining after removing rows with NA in effect columns.",
         call. = FALSE)

  data_file <- file.path(dir, "data_from_r.txt")
  utils::write.table(out_data, data_file, row.names = FALSE,
                      col.names = FALSE, quote = FALSE)

  # Call renumf90() with mapped positions
  renumf90(
    datafile = basename(data_file),
    traits = trait_positions,
    residual_variance = residual_variance,
    effects = effects,
    weights = weight_pos,
    snp_file = snp_file,
    inbreeding = inbreeding,
    ped_depth = ped_depth,
    options = options,
    dir = dir
  )
}


## Expand a variance specification to a proper n x n matrix.
## Handles: NULL (use default), scalar (diagonal), vector (diagonal),
## and matrix (pass through). Validates dimensions.
expand_variance <- function(x, n_traits, default_diag) {
  if (is.null(x)) {
    return(diag(default_diag, n_traits))
  }

  x <- as.matrix(x)

  # Scalar -> diagonal matrix
  if (length(x) == 1) {
    return(diag(as.numeric(x), n_traits))
  }

  # Vector -> diagonal matrix
  if (min(dim(x)) == 1 && max(dim(x)) == n_traits) {
    return(diag(as.numeric(x), n_traits))
  }

  # Matrix -> validate dimensions
  if (nrow(x) == n_traits && ncol(x) == n_traits) {
    return(x)
  }

  stop("Variance specification must be a scalar, vector of length ",
       n_traits, ", or a ", n_traits, "x", n_traits, " matrix.",
       call. = FALSE)
}
