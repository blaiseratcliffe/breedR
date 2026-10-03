# breedR examples

Worked examples beyond the [quick example in the README](../README.md#a-quick-example). Examples on `globulus` and `douglas` run as they are, since both data sets come with breedR. The others stand for your own data. Wherever an example uses an object it does not create, such as `dat`, `ped`, `renum_output` or `geno_matrix`, that object is yours, with the columns the example names. In the pedigree examples, `dat` has `self`, `sire` and `dam` in its first three columns. Names such as `genotypes.txt` or `pedigree.txt` are files you supply.

Several examples need programs beyond BLUPF90+: install them with `install_genomic_programs()` and `install_renumf90()`.

## Contents

- [Model fitting](#model-fitting)
  - [Spatial model with a grid search over rho](#spatial-model-with-a-grid-search-over-rho)
  - [Multi-trait models with different effects per trait](#multi-trait-models-with-different-effects-per-trait)
  - [Heterogeneous residual variances](#heterogeneous-residual-variances)
  - [Following and resuming a fit](#following-and-resuming-a-fit)
- [Heritability and genetic correlations](#heritability-and-genetic-correlations)
- [Genomic evaluation](#genomic-evaluation)
  - [Single-step GBLUP](#single-step-gblup)
  - [GWAS](#gwas)
  - [Predicting new animals](#predicting-new-animals)
- [Bayesian inference](#bayesian-inference)
- [Data preparation with RENUMF90](#data-preparation-with-renumf90)
- [Genotype QC and parentage](#genotype-qc-and-parentage)
- [Prediction validation](#prediction-validation)
- [SNP file I/O](#snp-file-io)
- [Reference tables](#reference-tables)

## Model fitting

### Spatial model with a grid search over rho

Without `rho`, an AR spatial model is fitted once for each pair on a grid of `rho` values, and the most likely fit is returned. `parallel` spreads the grid over several cores.

```r
library(breedR)
data(globulus)

res <- remlf90(
  fixed    = phe_X ~ gg,
  spatial  = list(model = 'AR',
                  coord = globulus[, c('x', 'y')]),
  data     = globulus,
  parallel = 4   # use 4 cores for the grid search
)

res$rho   # the grid, with the log-likelihood of each fit
```

### Multi-trait models with different effects per trait

When only the random effects differ between traits, fit the model with `remlf90()` and use `traits` to name the responses each random effect group applies to. Groups that `traits` does not list apply to every response. `traits` covers random formula terms, `genetic`, `spatial`, `pec` and generic groups (see `?remlf90`). With `genomic`, the genetic group must stay on every response: restricting it to some traits is refused.

In this example on the `douglas` trial, the block effect is fitted on circumference (`C13`) only, while the additive genetic effect covers both height (`H05`) and circumference:

```r
data(douglas)
dat <- droplevels(subset(douglas, site == "s1"))

res <- remlf90(
  fixed   = cbind(H05, C13) ~ orig,
  random  = ~ block,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, 1:3],
                 id = 'self'),
  traits  = list(block = "C13"),   # block on C13 only; genetic on both
  data    = dat
)

summary(res)            # reports "Absent random effects: block on H05."
res$var[["block", 1]]   # the H05 row and column are NA: not fitted
```

`traits` does not apply to fixed effects. When a fixed effect belongs to some traits only, as `age` does for height below, use `renumf90()` and give the effect position 0 for the traits it does not affect:

```r
# Height depends on herd + age; diameter depends on herd + site
renum <- renumf90(
  datafile = "data.txt",
  traits   = c(5, 6),    # height col 5, diameter col 6
  residual_variance = matrix(c(10, 3, 3, 5), 2, 2),
  effects  = list(
    list(pos = c(2, 2), type = "cross", form = "alpha"),  # herd: both
    list(pos = c(3, 0), type = "cov"),                     # age: height only
    list(pos = c(0, 4), type = "cross", form = "alpha"),  # site: diameter only
    list(pos = c(1, 1), type = "cross", form = "alpha",
         random = "animal", file = "pedigree.txt",
         covariances = matrix(c(5, 2, 2, 3), 2, 2))       # genetic: both
  )
)

# Fit with heritabilities and the genetic correlation
res <- remlf90_from_renum(renum, method = 'ai',
  progsf90.options = c(
    h2_formula(genetic_effect = 4, trait = 1),
    h2_formula(genetic_effect = 4, trait = 2),
    rg_formula(1, 2, effect = 4)
  )
)
```

[inst/doc/multi_trait_different_effects.md](../inst/doc/multi_trait_different_effects.md) walks through both routes in more detail.

### Heterogeneous residual variances

`hetres_options(covariate_cols = ...)` models the residual variance as a log-linear function of covariates, var(e) = exp(a0 + a1 x1 + ...). `covariate_cols` gives column positions in the data file that breedR writes for the backend. In this model its columns are `phe_X`, the intercept and `xc`, so `xc` is column 3. For another model, fit it once with `debug = TRUE` to see the file. This model needs AI-REML, the default.

```r
dat <- globulus
dat$xc <- dat$x / 100   # rescale so that a starting slope of 0.1 is sensible

res <- remlf90(
  fixed   = phe_X ~ xc,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, 1:3],
                 id = 'self'),
  data    = dat,
  progsf90.options = hetres_options(
    covariate_cols = 3,         # data file columns: phe_X, intercept, xc, ...
    initial = c(log(5), 0.1)    # a0, a1; a slope that starts at 0 never moves
  )
)

res$hetres   # a0 and a1, with standard errors
```

The fit has no `Residual` row in `res$var`, and no default heritability, which BLUPF90+ cannot compute with heterogeneous residuals. A residual variance per class of a factor is estimated by `gibbsf90(hetres_int = ...)`, not by `remlf90()`: see `?hetres_options` and the `Heterogeneous-variances` vignette.

### Following and resuming a fit

`progress_file` writes the backend's output to a file as it runs, so a long fit can be followed with `tail -F` or read with `reml_checkpoint()` from another R session. `cont = TRUE` resumes a fit from the last complete round in that file, for example after an interruption or a fit that stopped at `maxrounds`.

Both need a single local fit. An AR `rho` grid search (no `rho`, or a `rho` matrix) refuses them, so fix `rho` to one pair first. Remote fits (`breedR.bin = 'remote'` or `'submit'`) refuse them too. `cont = TRUE` takes the starting values from the file, so it refuses an explicit `var.ini` and the `initial` argument of `hetres_options()`. It also refuses `debug = TRUE`, which parses no results.

```r
fit_globulus <- function(...)
  remlf90(
    fixed   = phe_X ~ gg,
    genetic = list(model = 'add_animal',
                   pedigree = globulus[, 1:3],
                   id = 'self'),
    spatial = list(model = 'AR',
                   coord = globulus[, c('x', 'y')],
                   rho = c(.85, .8)),   # one pair: a grid is refused
    data = globulus,
    ...
  )

# Stop after 3 rounds, standing in for a fit that was interrupted
res <- fit_globulus(progress_file = "globulus.log",
                    progsf90.options = "maxrounds 3")
reml_checkpoint("globulus.log")   # variances at the last complete round

# Resume from there; the previous log is kept with a timestamp
res <- fit_globulus(progress_file = "globulus.log", cont = TRUE)
res$reml$resumed_from   # 3
```

## Heritability and genetic correlations

`remlf90()` adds a heritability to models with a genetic effect. The helpers build other functions of the variance components, with standard errors, for `progsf90.options`. Effects are numbered in model order: in this model `gg` is effect 1, the genetic effect 2 and the spatial effect 3.

```r
res <- remlf90(phe_X ~ gg,
  genetic = list(model = 'add_animal', pedigree = globulus[, 1:3], id = 'self'),
  spatial = list(model = 'AR', coord = globulus[, c('x', 'y')], rho = c(.85, .8)),
  data = globulus,
  progsf90.options = var_functions(n_traits = 1, genetic_effect = 2,
                                   other_random = 3)
)
# res$funvars has one column per function, with rows 'mean' (the estimate at
# the REML solution), 'sample mean' and 'sample sd' (mean and SE of the Monte
# Carlo draws). summary() prints all three.

# For three traits: heritabilities and all genetic correlations
var_functions(n_traits = 3, genetic_effect = 2, correlations = TRUE)
```

`h2_formula()`, `rg_formula()`, `vp_formula()` and `var_ratio_formula()` build single functions; the `Heritability` vignette explains the computation.

## Genomic evaluation

These are short versions. The `Genomic-selection` vignette, `vignette("Genomic-selection", package = "breedR")`, runs the whole pipeline from genotype preparation to prediction on example data, and [inst/doc/pregsf90_postgsf90_options.md](../inst/doc/pregsf90_postgsf90_options.md) lists the PREGSF90 and POSTGSF90 options.

### Single-step GBLUP

```r
res <- remlf90(
  fixed   = phe_X ~ gg,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, 1:3],
                 id = 'self'),
  data    = dat,
  genomic = list(
    snp_file = "genotypes.txt",
    map_file = "snp_map.txt",
    whichG   = 1,        # VanRaden (2008)
    tunedG   = 2,        # scale G to match A22
    minfreq  = 0.05      # MAF threshold
  )
)

# QC report from PREGSF90
res$genomic$n_snp_passed
res$genomic$excluded_snp
```

### GWAS

`postgsf90()` back-solves SNP effects from the genomic G-inverse, so fit the model with `save_ginverse = TRUE`. Without it, `postgsf90()` stops and asks you to refit.

```r
res <- remlf90(
  fixed   = phe_X ~ gg,
  genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),
  data    = dat,
  genomic = list(snp_file = "genotypes.txt",
                 map_file = "snp_map.txt",
                 save_ginverse = TRUE)     # required for a later GWAS
)

gwas <- postgsf90(res,
  windows_variance = 20,
  manhattan_plot   = TRUE
)

gwas$snp_sol      # SNP solutions and weights
gwas$manhattan    # Manhattan plot data
gwas$windows      # variance explained by genomic windows
```

### Predicting new animals

`predf90()` reads the SNP effects (`snp_pred`) that `postgsf90()` wrote. Pass it the directory they are in, `gwas$dir`. That is a temporary directory of the R session that fitted the model, and it is deleted when the session ends, so run `remlf90()`, `postgsf90()` and `predf90()` in one session. A saved `gwas` object reloaded later keeps the path but not the files.

```r
predictions <- predf90(
  snp_file   = "new_genotypes.txt",
  dir        = gwas$dir,   # where postgsf90() wrote snp_pred
  use_mu_hat = TRUE        # add the base so DGV are comparable to GEBV
)
# Returns: data.frame(id, call_rate, dgv)
```

Reliabilities (`acc = TRUE`) also need an `OPTION snp_var` file: see [Known limitations](../README.md#known-limitations).

## Bayesian inference

`gibbsf90()` takes `fixed`, `genetic`, `spatial` and the other model arguments of `remlf90()`, and estimates the variance components by Gibbs sampling. `postgibbsf90()` checks convergence.

```r
res <- gibbsf90(phe_X ~ gg,
  genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),
  data = dat, n_samples = 50000, burnin = 5000, thin = 10)

# Posterior means of the variance components
colMeans(res$samples)

# Convergence diagnostics
pg <- postgibbsf90(res, burnin = 5000, thin = 10)
pg$effective_size    # should be > 10
pg$geweke            # should be |x| < 2
pg$hpd               # 95% HPD intervals

# Threshold model for binary survival, coded 1 and 2 (0 = missing)
res_bin <- gibbsf90(survival ~ site,
  genetic = list(model = 'add_animal', pedigree = ped, id = 'self'),
  data = dat, n_samples = 100000, burnin = 20000, thin = 10,
  cat = 2)   # 2 = binary trait
```

`gibbsf90_from_renum()` runs the same sampler on `renumf90()` output.

## Data preparation with RENUMF90

RENUMF90 renumbers data files, checks the pedigree and sets up unknown parent groups. It suits data sets with alphanumeric ids. `renumf90()` works on column positions in a data file:

```r
renum <- renumf90(
  datafile = "raw_data.txt",
  traits   = c(5, 6),
  residual_variance = matrix(c(5, 2, 2, 4), 2, 2),
  effects  = list(
    list(pos = c(1, 1), type = "cross", form = "alpha"),
    list(pos = c(2, 2), type = "cross", form = "alpha",
         random = "animal",
         file = "pedigree.txt",
         covariances = matrix(c(10, 3, 3, 11), 2, 2))
  ),
  inbreeding = "pedigree"
)

# Fit the model on the RENUMF90 output
res <- remlf90_from_renum(renum, method = 'ai')
```

`renumf90_from_data()` takes an R data frame and column names instead, and writes the files itself:

```r
renum <- renumf90_from_data(
  data     = dat,
  traits   = "height",
  fixed    = c("site", "age"),
  random   = c("block"),
  pedigree = dat[, c("tree_id", "sire", "dam")],
  factors  = c("block"),    # treat a numeric block code as a factor
  genetic_variance  = 5,
  residual_variance = 10
)
res <- remlf90_from_renum(renum, method = 'ai')
```

`remlf90_from_renum()` returns a plain list of raw results, not a `remlf90` object, so `summary()` and the other methods do not apply to it.

## Genotype QC and parentage

`qcf90()` checks genotypes and pedigree on raw, not yet renumbered, data:

```r
qc <- qcf90(
  snp_file = "genotypes.txt",
  ped_file = "pedigree.txt",
  map_file = "snp_map.txt",
  maf = 0.05,
  check_parentage = TRUE,
  save_clean = TRUE
)
qc$log              # QC summary
qc$clean_snp        # path to the cleaned genotype file
qc$removed_animals  # animals excluded
```

`seekparentf90()` checks parent-offspring pairs against the SNPs and searches for the right parent when one conflicts:

```r
result <- seekparentf90(
  snp_file  = "genotypes.txt",
  ped_file  = "pedigree.txt",
  seek_sire = TRUE,      # search for the correct sire
  seek_dam  = TRUE,      # search for the correct dam
  yob       = TRUE       # year of birth in column 4 of the pedigree file
)
result$check       # pedigree with conflicting parents removed
result$assigned    # corrected pedigree after parent assignment
```

## Prediction validation

`validate_prediction()` runs the LR validation method (Legarra and Reverter, 2018): it fits the model on the whole data and on a copy where the validation animals' phenotypes for the focal trait are set to missing, and compares the two sets of predictions. `validation_ids` are RENUMF90's renumbered ids, not the original ones. The `pedigree` element of the `renumf90()` result maps one to the other: it is the pedigree of the first animal effect, with `animal` and `original_id` columns when RENUMF90 wrote its usual 10-column file. An id that is not found there would silently stay in the partial data, so the lookup below stops on one. The bundled `validationf90` cannot finish this yet: see [Known limitations](../README.md#known-limitations).

```r
# Original ids to renumbered ones
ped <- renum_output$pedigree
young_animal_ids <- ped$animal[match(as.character(young_original_ids),
                                     as.character(ped$original_id))]
stopifnot(!anyNA(young_animal_ids))   # every id must be in the pedigree

val <- validate_prediction(
  renum          = renum_output,     # from renumf90()
  validation_ids = young_animal_ids, # renumbered ids of the animals to validate
  effect         = 2,                # genetic effect number
  trait          = 1                 # focal trait to hold out (multi-trait)
)
val$statistics     # bias, dispersion, accuracy
```

## SNP file I/O

```r
# Write a genotype matrix in BLUPF90 format
write_snp_file(geno_matrix, animal_ids, "genotypes.txt")

# Fractional genotypes, for example from imputation
write_snp_file(imputed_geno, animal_ids, "genotypes.txt", fractional = TRUE)

# Read it back
snp_data <- read_snp_file("genotypes.txt")
snp_data$ids    # animal ids
snp_data$geno   # genotype matrix
```

## Reference tables

### BLUPF90 programs wrapped

| R function | BLUPF90 program | Purpose |
|---|---|---|
| `remlf90()` | BLUPF90+ | Variance component estimation (AI-REML / EM-REML) |
| `remlf90(genomic = ...)` | PREGSF90 + BLUPF90+ | Genomic QC, G matrix and ssGBLUP |
| `postgsf90()` | PostGSF90 | SNP effects, GWAS, Manhattan plots |
| `predf90()` | PREDF90 | Genomic prediction for new animals |
| `renumf90()` | RENUMF90 | Data renumbering, pedigree validation, UPGs |
| `remlf90_from_renum()` | BLUPF90+ | Model fitting from RENUMF90 output |
| `seekparentf90()` | SeekParentF90 | Parentage verification and discovery via SNP |
| `validationf90()` | ValidationF90 | LR prediction validation |
| `validate_prediction()` | BLUPF90+ + ValidationF90 | Full validation workflow |
| `gibbsf90()` | GIBBSF90+ | Bayesian variance estimation (Gibbs sampling) |
| `postgibbsf90()` | POSTGIBBSF90 | Convergence diagnostics for Gibbs samples |
| `gibbsf90_from_renum()` | GIBBSF90+ | Gibbs sampling from RENUMF90 output |
| `qcf90()` | QCF90 | Genotype and pedigree QC (before RENUMF90, raw ids) |

### Helper functions

| Function | Purpose |
|---|---|
| `h2_formula()` | Heritability option with SE |
| `rg_formula()` | Genetic correlation option with SE |
| `var_functions()` | All variance functions for a model |
| `vp_formula()` | Phenotypic variance option with SE |
| `var_ratio_formula()` | Variance proportion option with SE |
| `hetres_options()` | Heterogeneous residual variance options |
| `write_snp_file()` | Write a genotype matrix in BLUPF90 format |
| `read_snp_file()` | Read a BLUPF90 genotype file into R |
| `renumf90_from_data()` | Prepare data from R data frames, without column positions |