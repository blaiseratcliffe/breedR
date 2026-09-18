# Multi-Trait Models with Different Effects per Trait

Use `remlf90(traits = ...)` to include a random-effect group on selected
responses in a joint model. Fixed effects remain common to all responses.
For different fixed effects per trait (e.g., age affects height but not
diameter), use `renumf90()` + `remlf90_from_renum()`.

## Selecting random effects in the formula interface

This small example fits a replicate effect on height only:

```r
library(breedR)
a <- c(-3, -2, -1, 1, 2, 3)
dat <- data.frame(
  rep = factor(rep(seq_along(a), each = 4)),
  height = 10 + rep(a, each = 4) + rep(c(-1, 1, -1, 1), 6),
  diameter = 20 + rep(c(-2, -2, 2, 2), 6)
)
fit <- remlf90(
  cbind(height, diameter) ~ 1,
  random = ~ rep,
  traits = list(rep = "height"),
  var.ini = list(rep = diag(c(5, 0)), residuals = diag(c(1, 4))),
  data = dat
)
ranef(fit)$rep   # diameter coefficients are NA: this effect was not fitted
fitted(fit)      # diameter fitted values still include its fixed intercept
```

The `traits` list names random-effect groups and gives response-column names.
Unlisted groups apply to every response. `NULL`, an empty list, and explicitly
selecting every response give the usual model. Selection order does not change
response order; for example, `traits = list(rep = c("y3", "y1"))` selects the
first and third columns of `cbind(y1, y2, y3)`. A selection must retain at least
one response. Remove the group from the model to omit it everywhere.

Groups may be ordinary random formula terms, `genetic`, `spatial`, `pec`, or
named generic components. A correlated group shares one selection: direct and
competition genetic effects cannot have different selections. Generic groups
named `genetic` or `spatial` use their existing normalized names,
`generic_genetic` and `generic_spatial`.

Initial covariance matrices keep their full dimensions and existing order:
group member first, then response. The active principal block must be positive
definite. Inactive entries are ignored and their rows and columns are written
as zero; supply finite, symmetric numeric matrices without `NA`. For example,
both `diag(c(5, 0))` and `diag(c(5, 2))` initialize the fit above identically.
An active off-diagonal zero still constrains that covariance to zero. A zero
variance or an `NA` in `var.ini` does not select an absent effect.

Variance and coefficient results retain all response dimensions. `NA` marks
absent combinations; `summary()` identifies them. Their contribution to
`fitted()` is zero, and missing response observations stay missing in
`residuals()`. Model incidence matrices and observation order are unchanged.

Use the same `traits` specification when restarting with `cont = TRUE` and
`progress_file`. `reml_checkpoint(logfile, model = fit)` validates the active
covariance blocks and returns full matrices with zero rows for absent traits;
these can also be supplied in the existing `var.ini` locations with the same
selection. Pass the fitted model when reading such a checkpoint: zero entries
alone do not identify which effects were absent.

With `genomic`, restrictions on other random groups are supported, but the
genetic group must apply to all responses. Differing fixed effects and
member-specific selections within a correlated group require the lower-level
route below. That route returns raw results in a plain list, not a `remlf90`
object with `ranef()`, `fitted()`, and `summary()` methods.

## The Key Concept: `pos = c(col_trait1, col_trait2)`

Each effect in RENUMF90 gets a position vector — one column number per trait.
Use `0` to indicate an effect is absent for that trait.

```r
# Age covariate for trait 1 only (column 3 in data file):
list(pos = c(3, 0), type = "cov")

# Site factor for trait 2 only (column 4 in data file):
list(pos = c(0, 4), type = "cross", form = "alpha")

# Herd factor for both traits (column 2 in data file):
list(pos = c(2, 2), type = "cross", form = "alpha")
```

## Full Example: Height and Diameter with Different Models

Suppose you have forest trial data where:
- **Height** depends on herd + age + genetic
- **Diameter** depends on herd + site + genetic
- Both traits share the genetic effect (correlated breeding values)

### Step 1: Prepare Data

Your data file (`trial_data.txt`) should be space-separated with no header:

```
# animal_id  herd   age  site  height  diameter
tree001      herdA  5.2  s01   12.3    8.1
tree002      herdA  4.8  s01   11.7    7.5
tree003      herdB  6.1  s02   14.2    9.3
...
```

Your pedigree file (`pedigree.txt`):

```
# animal  sire  dam
tree001   sire1  dam1
tree002   sire1  dam2
tree003   sire2  dam1
sire1     0      0
sire2     0      0
dam1      0      0
dam2      0      0
```

If your data is in R as a data.frame:

```r
write.table(dat[, c("id", "herd", "age", "site", "height", "diameter")],
            "trial_data.txt", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
write.table(ped[, c("id", "sire", "dam")],
            "pedigree.txt", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
```

### Step 2: Run RENUMF90

```r
library(breedR)

renum <- renumf90(
  datafile = "trial_data.txt",
  traits = c(5, 6),       # height = col 5, diameter = col 6

  residual_variance = matrix(c(
    10, 3,                 # starting residual covariance
     3, 5
  ), 2, 2),

  effects = list(
    # Effect 1: Herd — both traits
    list(pos = c(2, 2), type = "cross", form = "alpha"),

    # Effect 2: Age — height (trait 1) only
    list(pos = c(3, 0), type = "cov"),

    # Effect 3: Site — diameter (trait 2) only
    list(pos = c(0, 4), type = "cross", form = "alpha"),

    # Effect 4: Animal — genetic effect, both traits
    list(pos = c(1, 1), type = "cross", form = "alpha",
         random = "animal",
         file = "pedigree.txt",
         covariances = matrix(c(
           5, 2,           # starting genetic covariance
           2, 3
         ), 2, 2))
  ),

  inbreeding = "pedigree"  # compute inbreeding coefficients
)
```

### Step 3: Check RENUMF90 Output

```r
# Generated parameter file (what BLUPF90+ will read)
cat(renum$par_content, sep = "\n")

# Should show something like:
# EFFECTS:
#  3  3   12 cross    <- herd: both traits
#  4  0    1 cov      <- age: trait 1 only
#  0  5    8 cross    <- site: trait 2 only
#  6  6  150 cross    <- animal: both traits

# Pedigree summary
cat("Animals in pedigree:", nrow(renum$pedigree), "\n")

# Inbreeding
head(renum$inbreeding)
```

### Step 4: Fit the Model

```r
res <- remlf90_from_renum(renum, method = 'ai')

# Solutions: breeding values, fixed effects
head(res$solutions)
```

### Step 5: Compute Genetic Parameters

Use the variance function helpers to get heritabilities and genetic
correlation with standard errors:

```r
# Heritability for height (trait 1) and diameter (trait 2)
# Effect 4 = genetic, effects 4 = only random effect
res <- remlf90_from_renum(renum, method = 'ai',
  progsf90.options = c(
    h2_formula(genetic_effect = 4, trait = 1),
    h2_formula(genetic_effect = 4, trait = 2),
    rg_formula(trait1 = 1, trait2 = 2, effect = 4)
  )
)
```

This will estimate:
- h2 for height = Va_height / (Va_height + Ve_height)
- h2 for diameter = Va_diameter / (Va_diameter + Ve_diameter)
- Genetic correlation = Cov_a(height, diameter) / sqrt(Va_height * Va_diameter)

All with standard errors via the Meyer & Houle (2013) sampling method.

## Three-Trait Example

For three traits with maximally different models:

```r
renum <- renumf90(
  datafile = "data.txt",
  traits = c(7, 8, 9),    # three traits in columns 7, 8, 9

  residual_variance = matrix(c(
    10, 3, 2,
     3, 8, 1,
     2, 1, 6
  ), 3, 3),

  effects = list(
    # Site: all three traits
    list(pos = c(2, 2, 2), type = "cross", form = "alpha"),

    # Block nested in site: traits 1 and 2 only
    list(pos = c(3, 3, 0), type = "cross", form = "alpha"),

    # Age covariate: trait 1 only
    list(pos = c(4, 0, 0), type = "cov"),

    # Density covariate: traits 2 and 3 only
    list(pos = c(0, 5, 5), type = "cov"),

    # Animal genetic: all three traits
    list(pos = c(1, 1, 1), type = "cross", form = "alpha",
         random = "animal",
         file = "pedigree.txt",
         covariances = matrix(c(
           5, 2, 1,
           2, 4, 1,
           1, 1, 3
         ), 3, 3))
  )
)

res <- remlf90_from_renum(renum, method = 'ai',
  progsf90.options = var_functions(
    n_traits = 3, genetic_effect = 5, correlations = TRUE
  )
)
```

## Adding Spatial Effects

For spatial models, include the coordinates as covariates and use the
AR or splines structure. With RENUMF90, spatial effects are specified
as generic random effects with a user-provided precision matrix:

```r
# Compute the AR precision matrix in R
ar_obj <- breedr_ar(
  coordinates = dat[, c("row", "col")],
  rho = c(0.85, 0.80)
)
precision_matrix <- vcov(ar_obj)

# Write the precision matrix in triplet format
triplet <- as.triplet(precision_matrix)
write.table(triplet, "ar_precision.txt",
            row.names = FALSE, col.names = FALSE)

# Add to the RENUMF90 effects list:
list(pos = c(6, 6), type = "cross", form = "numer",
     random = "diagonal",  # structure comes from user_file
     covariances = matrix(c(3, 1, 1, 2), 2, 2))
```

Note: integrating spatial effects through RENUMF90 requires manual
precision matrix construction. With common fixed effects, use `remlf90()`
and select the spatial group's responses with `traits = list(spatial = ...)`.

## When to Use Each Pathway

| Scenario | Use |
|---|---|
| Same fixed effects for all traits | `remlf90()` with `cbind(y1, y2) ~ ...` |
| Random-effect groups present on selected traits | `remlf90(traits = list(group = c("y1", "y3")))` |
| Different fixed effects per trait | `renumf90()` + `remlf90_from_renum()` |
| Alphanumeric IDs | `renumf90()` |
| Unknown parent groups | `renumf90()` |
| Simple single-trait model | `remlf90()` |
| Genomic ssGBLUP | `remlf90(genomic = ...)` |
