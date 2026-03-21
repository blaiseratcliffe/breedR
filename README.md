# breedR

### Statistical methods for genetic resources analysts — with genomic selection

breedR is an R package for fitting linear mixed models in breeding and quantitative genetics. It wraps the [BLUPF90](https://nce.ads.uga.edu/wiki/doku.php?id=start) family of programs to estimate variance components via REML, compute breeding values, and run genomic evaluations.

> **Note:** This fork is in active development and testing. The genomic features (ssGBLUP, GWAS, PREDF90, RENUMF90) are functional but have not been validated on production datasets. Use at your own risk and please report issues.

This fork extends the [original breedR](https://github.com/famuvie/breedR) with:

- **BLUPF90+ backend** — migrated from deprecated standalone binaries to the current unified program
- **Single-step genomic BLUP (ssGBLUP)** — PREGSF90 integration for combining pedigree and genomic information
- **Genome-wide association (ssGWAS)** — PostGSF90 wrapper for SNP effect extraction and Manhattan plots
- **Genomic prediction** — PREDF90 wrapper for predicting DGV in new animals
- **RENUMF90 data preparation** — full wrapper for data renumbering, pedigree validation, and unknown parent groups
- **Parentage verification** — SeekParentF90 wrapper for SNP-based paternity testing and parent discovery
- **Prediction validation** — ValidationF90 wrapper for LR validation method (Legarra & Reverter)
- **Variance function helpers** — `h2_formula()`, `rg_formula()`, `var_functions()` for heritability and genetic correlations with SEs
- **Heterogeneous residual variances** — `hetres_options()` for log-linear residual variance models
- **Bayesian inference** — GIBBSF90+ wrapper for Gibbs sampling, threshold/categorical traits, convergence diagnostics via POSTGIBBSF90
- **Genotype QC** — QCF90 wrapper for standalone quality control on raw (non-renumbered) data
- **SNP file I/O** — `write_snp_file()` / `read_snp_file()` for integer and fractional genotype formats
- **Convenience data preparation** — `renumf90_from_data()` takes R data.frames and column names instead of file positions
- **Performance improvements** — vectorized spatial indexing, parallel AR grid search, resilient error handling

## Installation

```r
# Install from this fork
devtools::install_github('blaiseratcliffe/breedR')

# Or install from source
git clone https://github.com/blaiseratcliffe/breedR.git
R CMD INSTALL breedR
```

BLUPF90+ binaries are downloaded automatically from UGA at install time. To install additional programs:

```r
library(breedR)
install_genomic_programs()   # preGSf90, postGSf90, predf90, validationf90,
                             # predictf90, seekparentf90, gibbsf90+,
                             # postgibbsf90, qcf90
install_renumf90()           # renumf90
```

## Quick Start

### Basic animal model

```r
library(breedR)
data(globulus)

res <- remlf90(
  fixed   = phe_X ~ gg,
  genetic = list(model = 'add_animal',
                 pedigree = globulus[, 1:3],
                 id = 'self'),
  spatial = list(model = 'AR',
                 coord = globulus[, c('x', 'y')],
                 rho = c(.85, .8)),
  data = globulus
)

summary(res)
```

### Spatial model with automatic rho selection

```r
# Grid search over rho values (parallelized on multi-core systems)
res <- remlf90(
  fixed   = phe_X ~ gg,
  spatial = list(model = 'AR',
                 coord = globulus[, c('x', 'y')]),
  data    = globulus,
  parallel = 4   # use 4 cores for grid search
)

res$rho   # evaluation grid with log-likelihoods
```

### Genomic evaluation (ssGBLUP)

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

# QC report
res$genomic$n_snp_passed
res$genomic$excluded_snp
```

### GWAS

```r
gwas <- postgsf90(res,
  windows_variance = 20,
  manhattan_plot    = TRUE
)

gwas$snp_sol      # SNP solutions and weights
gwas$manhattan    # Manhattan plot data
gwas$windows      # variance explained by genomic windows
```

### Bayesian inference (Gibbs sampling)

```r
# Linear model with Gibbs sampling
res <- gibbsf90(phe_X ~ gg,
  genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),
  data = dat, n_samples = 50000, burnin = 5000, thin = 10)

# Posterior means for variance components
colMeans(res$samples)

# Convergence diagnostics
diag <- postgibbsf90(res, burnin = 5000, thin = 10)
diag$effective_size    # should be > 10
diag$geweke            # should be |x| < 2
diag$hpd               # 95% HPD intervals

# Threshold model for binary survival
res_bin <- gibbsf90(survival ~ site,
  genetic = list(model = 'add_animal', pedigree = ped, id = 'self'),
  data = dat, n_samples = 100000, burnin = 20000, thin = 10,
  cat = 2)   # 2 = binary trait
```

### Predict new animals

```r
predictions <- predf90(
  snp_file     = "new_genotypes.txt",
  use_mu_hat   = TRUE,
  acc          = TRUE
)
# Returns: data.frame(id, call_rate, dgv, reliability)
```

### RENUMF90 data preparation

```r
# For datasets with alphanumeric IDs, unknown parent groups, etc.
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

# Fit model using RENUMF90 output
res <- remlf90_from_renum(renum, method = 'ai')
```

### Easy data preparation (from R data.frames)

```r
# No need to think about column positions — just use column names
renum <- renumf90_from_data(
  data     = dat,
  traits   = "height",
  fixed    = c("site", "age"),
  random   = c("block"),
  pedigree = dat[, c("tree_id", "sire", "dam")],
  factors  = c("block"),    # treat numeric block as factor
  genetic_variance  = 5,
  residual_variance = 10
)
res <- remlf90_from_renum(renum, method = 'ai')
```

### Genotype quality control

```r
# QCF90 runs on raw (non-renumbered) data before RENUMF90
qc <- qcf90(
  snp_file = "genotypes.txt",
  ped_file = "pedigree.txt",
  map_file = "snp_map.txt",
  maf = 0.05,
  check_parentage = TRUE,
  save_clean = TRUE
)
qc$log              # QC summary
qc$clean_snp        # path to cleaned genotype file
qc$removed_animals  # animals excluded
```

### Multi-trait with different effects per trait

```r
# Height depends on herd + age; diameter depends on herd + site
# Use renumf90() with per-trait positions (0 = absent for that trait)
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

# Fit with heritabilities + genetic correlation
res <- remlf90_from_renum(renum, method = 'ai',
  progsf90.options = c(
    h2_formula(genetic_effect = 4, trait = 1),
    h2_formula(genetic_effect = 4, trait = 2),
    rg_formula(1, 2, effect = 4)
  )
)
```

See `inst/doc/multi_trait_different_effects.md` for the full workflow.

### Heritability and genetic correlations

```r
# Automatic heritability + SE for a model with genetic (effect 2) + spatial (effect 3)
res <- remlf90(phe_X ~ gg,
  genetic = list(model = 'add_animal', pedigree = dat[, 1:3], id = 'self'),
  spatial = list(model = 'AR', coord = dat[, c('x', 'y')], rho = c(.85, .8)),
  data = dat,
  progsf90.options = var_functions(n_traits = 1, genetic_effect = 2,
                                    other_random = 3)
)
# h2 and spatial proportion with SEs in res$funvars

# For multi-trait: heritabilities + all genetic correlations
var_functions(n_traits = 3, genetic_effect = 2, correlations = TRUE)
```

### Parentage verification

```r
result <- seekparentf90(
  snp_file = "genotypes.txt",
  ped_file = "pedigree.txt",
  seek_sire = TRUE,      # search for correct sire
  seek_dam  = TRUE,      # search for correct dam
  yob       = TRUE       # use year of birth from pedigree col 4
)
result$check       # Match/No-Match for each parent-offspring pair
result$assigned    # corrected pedigree after parent assignment
```

### Prediction validation

```r
# Full LR validation workflow (Legarra & Reverter 2018)
val <- validate_prediction(
  renum            = renum_output,     # from renumf90()
  validation_ids   = young_animal_ids, # animals to validate
  effect           = 2                 # genetic effect number
)
val$statistics     # bias, dispersion, accuracy
```

### SNP file formatting

```r
# Write genotype matrix in BLUPF90 format
write_snp_file(geno_matrix, animal_ids, "genotypes.txt")

# Fractional genotypes (from imputation)
write_snp_file(imputed_geno, animal_ids, "genotypes.txt", fractional = TRUE)

# Read back
snp_data <- read_snp_file("genotypes.txt")
snp_data$ids    # animal IDs
snp_data$geno   # genotype matrix
```

## Supported Models

| Model | Description |
|---|---|
| Animal model | Additive genetic effects via pedigree |
| Competition model | Direct + competition genetic effects with neighbour structure |
| AR(1)xAR(1) spatial | Autoregressive spatial correlation (row x column) |
| B-splines spatial | Bidimensional penalized splines |
| Blocks | Block/group random effects |
| Generic | User-supplied incidence + covariance/precision matrices |
| Multi-trait | Multiple correlated traits |
| Threshold/categorical | Binary (survival) and ordinal traits via Gibbs sampling |
| ssGBLUP | Single-step genomic BLUP combining pedigree + genomic data |

## BLUPF90 Programs Wrapped

| R Function | BLUPF90 Program | Purpose |
|---|---|---|
| `remlf90()` | BLUPF90+ | Variance component estimation (AI-REML / EM-REML) |
| `remlf90(genomic=...)` | PREGSF90 + BLUPF90+ | Genomic QC + G matrix + ssGBLUP |
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
| `qcf90()` | QCF90 | Genotype/pedigree QC (pre-RENUMF90, raw IDs) |

## Helper Functions

| Function | Purpose |
|---|---|
| `h2_formula()` | Generate heritability OPTION with SE |
| `rg_formula()` | Generate genetic correlation OPTION with SE |
| `var_functions()` | Generate all variance functions for a model |
| `vp_formula()` | Generate phenotypic variance OPTION with SE |
| `var_ratio_formula()` | Generate variance proportion OPTION with SE |
| `hetres_options()` | Generate heterogeneous residual variance OPTIONs |
| `write_snp_file()` | Write genotype matrix in BLUPF90 format |
| `read_snp_file()` | Read BLUPF90 genotype file into R |
| `write_xref_file()` | Write cross-reference ID file |
| `renumf90_from_data()` | Prepare data from R data.frames (no column positions needed) |
| `build_gibbs_options()` | Generate GIBBSF90+ OPTION strings |

## Requirements

- R >= 3.1.2
- Key dependencies: `Matrix`, `sp`, `ggplot2`, `pedigree`, `pedigreemm`
- BLUPF90+ binary (downloaded automatically from [UGA](https://nce.ads.uga.edu/html/projects/programs/))

## Development

```bash
# Run tests
Rscript -e "testthat::test_dir('tests/testthat')"

# Regenerate documentation
Rscript -e "roxygen2::roxygenise()"

# Build package
R CMD build . --no-build-vignettes
R CMD check breedR_*.tar.gz --no-vignettes
```

## Credits

breedR was originally developed by [Facundo Munoz](https://github.com/famuvie) as part of the Trees4Future and ProCoGen projects. The BLUPF90 programs are developed by [Ignacy Misztal's group](https://nce.ads.uga.edu/) at the University of Georgia.

## License

GPL-3
