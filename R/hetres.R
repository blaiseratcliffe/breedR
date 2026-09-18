### Heterogeneous residual variance helpers ###


#' Generate BLUPF90+ options for heterogeneous residual variances
#'
#' Creates the OPTION strings needed to model heterogeneous residual
#' variances. BLUPF90+ supports two approaches: class-based (different
#' variance per group) and polynomial (variance as a function of covariates).
#' Only the polynomial approach is estimated by \code{\link{remlf90}}.
#'
#' @param group_col integer. Column position in the data file containing
#'   the group indicator (for class-based heterogeneity). Used with
#'   \code{n_groups}.
#' @param n_groups integer. Number of groups/classes for heterogeneous
#'   variances. Required when \code{group_col} is specified.
#' @param covariate_cols integer vector. Column positions of covariates
#'   for polynomial residual heterogeneity (for covariate-based). The
#'   residual variance is modeled as exp(a0 + a1*X1 + a2*X2 + ...).
#' @param initial numeric vector. Initial values for the polynomial
#'   coefficients. For class-based: not needed (read from a file).
#'   For covariate-based: the intercept a0 followed by regression
#'   coefficients a1, a2, etc. Start the slopes at small non-zero values,
#'   e.g. \code{c(log(residual_var), rep(0.01, n_covariates))}. A coefficient
#'   that starts at exactly 0 is never updated, and BLUPF90+ (version 2.73)
#'   then crashes.
#' @param var_file character or NULL. Path to a file containing the
#'   residual (co)variances for each class. Required with \code{group_col}.
#'
#' @return Character vector of OPTION strings to pass to
#'   \code{progsf90.options}. The covariate-based ones are for
#'   \code{\link{remlf90}}.
#'
#' @details
#' \strong{Note:} With \code{\link{remlf90}}, use the covariate-based approach
#' (\code{covariate_cols}). BLUPF90+ reads the class-based options
#' (\code{group_col}) under AI-REML and EM-REML alike, echoes the class
#' variances, and then fits a homoscedastic model, so \code{remlf90},
#' \code{\link{remlf90_from_renum}} and \code{\link{validate_prediction}}
#' refuse them. To estimate a residual variance per class, use
#' \code{\link{gibbsf90}} with its \code{hetres_int} argument. The
#' \code{se_covar_function} for
#' heritability is not available with heterogeneous residuals, so
#' \code{\link{remlf90}} does not add its default heritability to such a fit.
#' The covariate-based model needs AI-REML: BLUPF90+ refuses it under EM-REML.
#' \code{remlf90} returns the estimated coefficients, with standard errors
#' from the inverse AI matrix, in the \code{hetres} element of the fit, and
#' the \code{var} element then has no \code{Residual} row.
#'
#' Two types of heterogeneous residual models are supported:
#'
#' \strong{Class-based} (\code{hetres_int}): Different residual variance
#' for each level of a grouping factor (e.g., site, herd, age class).
#' Uses \code{group_col} and \code{n_groups}. Estimated by Gibbs sampling
#' only (see above).
#'
#' \strong{Covariate-based} (\code{hetres_pos/hetres_pol}): Residual
#' variance modeled as a smooth function of one or more covariates via
#' log-linear regression: \eqn{\sigma^2_e = \exp(a_0 + a_1 X_1 + ...)}.
#' Uses \code{covariate_cols} and \code{initial}.
#'
#' For multi-trait models, covariate positions and initial values should
#' be specified trait-first: e.g., for 2 traits with 1 covariate each,
#' \code{covariate_cols = c(10, 10)} and
#' \code{initial = c(4.0, 4.0, 0.1, 0.1)}.
#'
#' @examples
#' \dontrun{
#' # Class-based: different residual variance per site (column 5, 10 sites),
#' # estimated by Gibbs sampling. remlf90() refuses class-based options.
#' gibbsf90(phe ~ fixed,
#'          genetic = list(...),
#'          data = dat,
#'          hetres_int = list(col = 5, n = 10))
#'
#' # Covariate-based: residual variance as function of age (column 8)
#' remlf90(phe ~ fixed,
#'         data = dat,
#'         progsf90.options = hetres_options(
#'           covariate_cols = 8,
#'           initial = c(log(10), 0.1)   # intercept + slope
#'         ))
#' }
#' @export
hetres_options <- function(group_col = NULL,
                            n_groups = NULL,
                            covariate_cols = NULL,
                            initial = NULL,
                            var_file = NULL) {

  opts <- character(0)

  if (!is.null(group_col) && !is.null(covariate_cols))
    stop("Specify either 'group_col' (class-based) or 'covariate_cols' ",
         "(covariate-based), not both.", call. = FALSE)

  if (!is.null(group_col)) {
    # Class-based heterogeneous residuals. BLUPF90+ reads the initial per-class
    # residual variances from the file named by hetres_var. remlf90() does not
    # synthesize this file, so it must be supplied; otherwise the option is
    # written but the backend has no initial variances and silently fails.
    if (is.null(n_groups))
      stop("'n_groups' is required with 'group_col'.", call. = FALSE)
    if (is.null(var_file))
      stop("'var_file' is required with 'group_col': class-based heterogeneous ",
           "residuals need a file of initial per-class variances.",
           call. = FALSE)
    opts <- c(opts, paste("hetres_int", group_col, n_groups))
    opts <- c(opts, paste("hetres_var", var_file))

  } else if (!is.null(covariate_cols)) {
    # Covariate-based (polynomial) heterogeneous residuals
    opts <- c(opts, paste("hetres_pos",
                           paste(covariate_cols, collapse = " ")))
    if (!is.null(initial))
      opts <- c(opts, paste("hetres_pol",
                             paste(initial, collapse = " ")))
  } else {
    stop("Specify either 'group_col' or 'covariate_cols'.", call. = FALSE)
  }

  return(opts)
}


## Class-based heterogeneous residuals (OPTION hetres_int) are not a REML
## option of BLUPF90+. Under 'method VCE', which every breedR REML entry point
## writes, version 2.73 reads the option and its variance file, echoes the class
## variances as "Fixed R", and then fits and solves a homoscedastic model: the
## estimates, -2logL and solutions match a fit without it, under AI and EM-REML
## alike (#39). Refuse it rather than return that fit as if it were requested.
## `opts` may hold progsf90.options (no OPTION prefix) or parameter file lines.
refuse_hetres_int <- function(opts) {
  if (any(grepl('^\\s*(OPTION\\s+)?hetres_int\\b', opts)))
    stop("Class-based heterogeneous residual variances (OPTION hetres_int) ",
         "are not estimated by BLUPF90+ REML, which would fit a homoscedastic ",
         "model instead.\n",
         " Use gibbsf90(hetres_int = list(col = , n = )) to estimate a ",
         "residual variance per class, or hetres_options(covariate_cols = ...) ",
         "to model the residual variance as a log-linear function of ",
         "covariates.", call. = FALSE)
  invisible(NULL)
}
