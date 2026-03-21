### Heterogeneous residual variance helpers ###


#' Generate BLUPF90+ options for heterogeneous residual variances
#'
#' Creates the OPTION strings needed to model heterogeneous residual
#' variances in \code{\link{remlf90}}. BLUPF90+ supports two approaches:
#' class-based (different variance per group) and polynomial (variance
#' as a function of covariates).
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
#'   coefficients a1, a2, etc. A reasonable starting point is
#'   \code{c(log(residual_var), rep(0, n_covariates))}.
#' @param var_file character or NULL. Path to a file containing initial
#'   residual (co)variances for each class (used with \code{group_col}).
#'   If NULL, BLUPF90+ estimates them from the data.
#'
#' @return Character vector of OPTION strings to pass to
#'   \code{progsf90.options} in \code{\link{remlf90}}.
#'
#' @details
#' \strong{Note:} When using heterogeneous residual variances with AIREMLF90
#' (via BLUPF90+), use the covariate-based approach (\code{covariate_cols}).
#' The class-based approach (\code{group_col}) is designed for GIBBS3F90 and
#' requires an initial variance file. The \code{se_covar_function} for
#' heritability is not available with heterogeneous residuals.
#'
#' Two types of heterogeneous residual models are supported:
#'
#' \strong{Class-based} (\code{hetres_int}): Different residual variance
#' for each level of a grouping factor (e.g., site, herd, age class).
#' Uses \code{group_col} and \code{n_groups}.
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
#' # Class-based: different residual variance per site (column 5, 10 sites)
#' remlf90(phe ~ fixed,
#'         genetic = list(...),
#'         data = dat,
#'         progsf90.options = hetres_options(group_col = 5, n_groups = 10))
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
    # Class-based heterogeneous residuals.
    # BLUPF90+ reads initial per-class variances from a file called 'hetres'.
    # If no var_file is provided, we create one with equal initial variances.
    if (is.null(n_groups))
      stop("'n_groups' is required with 'group_col'.", call. = FALSE)
    opts <- c(opts, paste("hetres_int", group_col, n_groups))
    if (!is.null(var_file))
      opts <- c(opts, paste("hetres_var", var_file))
    # Store n_groups as an attribute so remlf90() can create the hetres file
    attr(opts, "hetres_n_groups") <- n_groups

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
