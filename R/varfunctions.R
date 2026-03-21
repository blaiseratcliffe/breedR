### Variance function helpers for BLUPF90+ ###
### Generates OPTION se_covar_function lines from model structure ###


#' Generate heritability formula for BLUPF90+
#'
#' Creates the \code{OPTION se_covar_function} string for heritability
#' estimation with standard errors. Handles single-trait, multi-trait,
#' and maternal effect models.
#'
#' @param genetic_effect integer. Effect number for the additive genetic
#'   effect in the BLUPF90+ parameter file (default 2 — typical for
#'   models with one fixed effect + one genetic effect).
#' @param trait integer. Trait number (default 1). For multi-trait models,
#'   specify which trait's heritability to compute.
#' @param maternal_effect integer or NULL. Effect number for the maternal
#'   genetic effect. If provided, computes direct heritability excluding
#'   maternal variance.
#' @param other_random integer vector or NULL. Effect numbers of other
#'   random effects (e.g., permanent environment, spatial) to include
#'   in the denominator.
#' @param label character. Label for the heritability estimate
#'   (default "H2" or "H2_t<trait>").
#' @return Character string for \code{progsf90.options}.
#'
#' @examples
#' \dontrun{
#' # Simple animal model: h2 = Va / (Va + Ve)
#' remlf90(..., progsf90.options = h2_formula())
#'
#' # With spatial effect (effect 3): h2 = Va / (Va + Vs + Ve)
#' remlf90(..., progsf90.options = h2_formula(other_random = 3))
#'
#' # Multi-trait, trait 2
#' remlf90(..., progsf90.options = h2_formula(trait = 2))
#'
#' # Maternal model: direct h2
#' remlf90(..., progsf90.options = h2_formula(maternal_effect = 3))
#' }
#' @export
h2_formula <- function(genetic_effect = 2L,
                        trait = 1L,
                        maternal_effect = NULL,
                        other_random = NULL,
                        label = NULL) {

  t <- as.integer(trait)
  g <- as.integer(genetic_effect)

  # Numerator: additive genetic variance for this trait
  numerator <- paste0("G_", g, "_", g, "_", t, "_", t)

  # Denominator: all variance components
  denom_parts <- numerator

  # Maternal genetic effect
  if (!is.null(maternal_effect)) {
    m <- as.integer(maternal_effect)
    # Direct-maternal covariance (counted once)
    denom_parts <- c(denom_parts,
                      paste0("G_", g, "_", m, "_", t, "_", t))
    # Maternal genetic variance
    denom_parts <- c(denom_parts,
                      paste0("G_", m, "_", m, "_", t, "_", t))
  }

  # Other random effects
  if (!is.null(other_random)) {
    for (r in as.integer(other_random)) {
      denom_parts <- c(denom_parts,
                        paste0("G_", r, "_", r, "_", t, "_", t))
    }
  }

  # Residual
  denom_parts <- c(denom_parts, paste0("R_", t, "_", t))

  denominator <- paste(denom_parts, collapse = "+")
  formula <- paste0(numerator, "/(", denominator, ")")

  if (is.null(label)) {
    label <- if (t == 1L) "H2" else paste0("H2_t", t)
  }

  paste("se_covar_function", label, formula)
}


#' Generate genetic correlation formula for BLUPF90+
#'
#' Creates the \code{OPTION se_covar_function} string for genetic
#' correlation between two traits with standard errors.
#'
#' @param trait1 integer. First trait number.
#' @param trait2 integer. Second trait number.
#' @param effect integer. Effect number for the genetic effect (default 2).
#' @param label character or NULL. Label for the correlation estimate.
#' @return Character string for \code{progsf90.options}.
#'
#' @examples
#' \dontrun{
#' # Genetic correlation between traits 1 and 2
#' remlf90(..., progsf90.options = rg_formula(1, 2))
#'
#' # Between traits 1 and 3 for effect 3
#' remlf90(..., progsf90.options = rg_formula(1, 3, effect = 3))
#' }
#' @export
rg_formula <- function(trait1, trait2, effect = 2L, label = NULL) {

  t1 <- as.integer(trait1)
  t2 <- as.integer(trait2)
  e <- as.integer(effect)

  # rg = cov(t1,t2) / sqrt(var(t1) * var(t2))
  covariance <- paste0("G_", e, "_", e, "_", t1, "_", t2)
  var1 <- paste0("G_", e, "_", e, "_", t1, "_", t1)
  var2 <- paste0("G_", e, "_", e, "_", t2, "_", t2)

  formula <- paste0(covariance, "/(", var1, "*", var2, ")**0.5")

  if (is.null(label))
    label <- paste0("rg_", t1, t2)

  paste("se_covar_function", label, formula)
}


#' Generate phenotypic variance formula for BLUPF90+
#'
#' Creates the \code{OPTION se_covar_function} string for total
#' phenotypic variance with standard errors.
#'
#' @param trait integer. Trait number (default 1).
#' @param random_effects integer vector. All random effect numbers
#'   in the model.
#' @param label character or NULL.
#' @return Character string for \code{progsf90.options}.
#' @export
vp_formula <- function(trait = 1L, random_effects = 2L, label = NULL) {

  t <- as.integer(trait)

  parts <- character(0)
  for (r in as.integer(random_effects)) {
    parts <- c(parts, paste0("G_", r, "_", r, "_", t, "_", t))
  }
  parts <- c(parts, paste0("R_", t, "_", t))

  formula <- paste(parts, collapse = "+")

  if (is.null(label))
    label <- if (t == 1L) "Vp" else paste0("Vp_t", t)

  paste("se_covar_function", label, formula)
}


#' Generate variance ratio formula for BLUPF90+
#'
#' Creates the \code{OPTION se_covar_function} string for the ratio
#' of any random effect variance to total phenotypic variance.
#'
#' @param effect integer. Effect number whose proportion to compute.
#' @param trait integer. Trait number (default 1).
#' @param all_random integer vector. All random effect numbers in the model.
#' @param label character or NULL.
#' @return Character string for \code{progsf90.options}.
#'
#' @examples
#' \dontrun{
#' # Proportion of spatial variance: Vs / Vp
#' remlf90(..., progsf90.options = var_ratio_formula(
#'   effect = 3, all_random = c(2, 3)))
#' }
#' @export
var_ratio_formula <- function(effect, trait = 1L, all_random = 2L,
                               label = NULL) {

  t <- as.integer(trait)
  e <- as.integer(effect)

  numerator <- paste0("G_", e, "_", e, "_", t, "_", t)

  denom_parts <- character(0)
  for (r in as.integer(all_random)) {
    denom_parts <- c(denom_parts, paste0("G_", r, "_", r, "_", t, "_", t))
  }
  denom_parts <- c(denom_parts, paste0("R_", t, "_", t))

  formula <- paste0(numerator, "/(", paste(denom_parts, collapse = "+"), ")")

  if (is.null(label))
    label <- paste0("prop_e", e)

  paste("se_covar_function", label, formula)
}


#' Generate multiple variance function formulas at once
#'
#' Convenience function that generates heritability, genetic correlations,
#' and variance proportions for a full model in one call.
#'
#' @param n_traits integer. Number of traits.
#' @param genetic_effect integer. Genetic effect number (default 2).
#' @param other_random integer vector or NULL. Other random effect numbers.
#' @param maternal_effect integer or NULL. Maternal effect number.
#' @param correlations logical. Include genetic correlations for
#'   multi-trait models (default TRUE).
#' @return Character vector of OPTION strings for \code{progsf90.options}.
#'
#' @examples
#' \dontrun{
#' # Single trait with genetic (2) and spatial (3) effects
#' remlf90(..., progsf90.options = var_functions(
#'   n_traits = 1, genetic_effect = 2, other_random = 3))
#'
#' # Two-trait model: heritabilities + genetic correlation
#' remlf90(..., progsf90.options = var_functions(n_traits = 2))
#' }
#' @export
var_functions <- function(n_traits = 1L,
                           genetic_effect = 2L,
                           other_random = NULL,
                           maternal_effect = NULL,
                           correlations = TRUE) {

  opts <- character(0)

  # Heritability for each trait
  for (t in seq_len(n_traits)) {
    opts <- c(opts, h2_formula(
      genetic_effect = genetic_effect,
      trait = t,
      maternal_effect = maternal_effect,
      other_random = other_random
    ))
  }

  # Genetic correlations between all trait pairs
  if (correlations && n_traits > 1) {
    for (t1 in seq_len(n_traits - 1)) {
      for (t2 in (t1 + 1):n_traits) {
        opts <- c(opts, rg_formula(t1, t2, effect = genetic_effect))
      }
    }
  }

  # Variance proportions for other random effects
  if (!is.null(other_random)) {
    all_random <- c(genetic_effect, other_random)
    for (r in other_random) {
      for (t in seq_len(n_traits)) {
        opts <- c(opts, var_ratio_formula(
          effect = r, trait = t, all_random = all_random
        ))
      }
    }
  }

  return(opts)
}
