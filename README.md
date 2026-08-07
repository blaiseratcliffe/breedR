# breedR

### Statistical methods for genetic resources analysts — with genomic selection

breedR is an R package for fitting linear mixed models in breeding and quantitative genetics. It wraps the [BLUPF90](https://nce.ads.uga.edu/wiki/doku.php?id=start) family of programs to estimate variance components via REML, compute breeding values, and run genomic evaluations.

> **Note:** This fork is in active development. The genomic paths (ssGBLUP, GWAS, genomic prediction, RENUMF90, parentage verification, and QC) now run end to end and are verified on small test datasets, but they have **not** been validated on production data. Use at your own risk and please report issues. See [Known limitations](#known-limitations) for two features that need a newer binary or extra setup.

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

## Known limitations

- **`validate_prediction()` / `validationf90`** — the breedR wrapper is verified, but the bundled
  `validationf90` v1.01 aborts on valid input (an end-of-file read error inside the program). LR
  validation therefore needs a newer `validationf90` build; this is a binary-side limitation, not a
  wrapper defect. The whole/partial model fits and the partial-data construction work correctly.
- **`predf90(acc = TRUE)`** — reliabilities require an `OPTION snp_var` file from the `postgsf90()`
  run, which is not wired up yet. Direct genomic values (`acc = FALSE`, the default) work.

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

> For a full walkthrough of the genomic pipeline (genotype prep → ssGBLUP → GWAS
> → prediction), see the *Genomic selection* vignette:
> `vignette("Genomic-selection", package = "breedR")`.

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

`postgsf90()` back-solves SNP effects from the genomic G-inverse, so the model must be fitted with
`save_ginverse = TRUE` (otherwise `postgsf90()` stops and asks you to refit):

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

`predf90()` reads the SNP effects (`snp_pred`) written by `postgsf90()`, so run it in the same R
session as the GWAS step above:

```r
# ... after gwas <- postgsf90(res, ...) in the same session ...
predictions <- predf90(
  snp_file   = "new_genotypes.txt",
  use_mu_hat = TRUE       # add the base so DGV are comparable to GEBV
)
# Returns: data.frame(id, call_rate, dgv)
```

Reliabilities (`acc = TRUE`) additionally require an `OPTION snp_var` file — see
[Known limitations](#known-limitations).

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

> Requires a `validationf90` build that runs to completion — see
> [Known limitations](#known-limitations).

```r
# Full LR validation workflow (Legarra & Reverter 2018)
val <- validate_prediction(
  renum            = renum_output,     # from renumf90()
  validation_ids   = young_animal_ids, # animals to validate
  effect           = 2,                # genetic effect number
  trait            = 1                 # focal trait to hold out (multi-trait)
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

## Comparison with other mixed-model tools

An at-a-glance comparison with eight widely used alternatives for
quantitative-genetics mixed models. Capabilities of all nine evolve — verify
the details against the current documentation for your use case.

| | **breedR** | **ASReml-R** | **sommer** | **lme4** | **gremlin** | **SpATS** | **BGLR** | **MCMCglmm** | **hibayes** |
|---|---|---|---|---|---|---|---|---|---|
| Engine | BLUPF90 (Fortran) | proprietary (VSNi) | R / C++ (Armadillo) | R / C++ (Eigen) | R / C++ (Matrix) | R (SAP algorithm) | R / C (Gibbs sampler) | R / C (MCMC) | R / C++ (MCMC) |
| License | GPL-3 (free) | commercial | GPL (free) | GPL (≥ 2) | GPL-3 (free) | GPL-2/3 (free) | GPL-3 (free) | GPL (≥ 2) | GPL-3 (free) |
| REML | Yes (AI / EM) | Yes | Yes | Yes (the default) | Yes (AI) | Yes (SAP) | — (Bayesian only) | — (Bayesian only) | — (Bayesian only) |
| Bayesian / MCMC | Yes (GIBBSF90+) | — | — | — | — | — | Yes (the entire package) | Yes (the entire package) | Yes (the entire package) |
| Pedigree / **A** matrix | Yes (built-in) | Yes | Yes | — (needs `pedigreemm`) | Yes (generalised inverse) | — | via RKHS kernel | Yes (`pedigree=`) | Yes (`ssbrm()`) |
| Spatial (AR, splines) | Yes | Yes | Yes | — | via user-supplied inverse | Yes (2D P-splines) | via user-supplied kernel | via `ginverse=` | — |
| Multi-trait | Yes | Yes | Yes | — | — (univariate only) | — | Yes (`Multitrait()`) | Yes (a core strength) | — |
| Single-step GBLUP | Yes (dedicated pipeline) | via user-supplied H | via user-supplied H | — | via user-supplied H⁻¹ | — | via user-supplied H kernel | via `ginverse=` (H⁻¹) | Yes (`ssbrm()`) |
| GWAS / SNP effects | Yes (PostGSF90) | via SNP models | Yes | — | — | — | Yes (BayesB/C + `BFDR()`) | — | Yes (PIP / WPPA) |
| Competition / IGE | Yes (built-in) | via custom structures | — | — | — | — | — | — | — |
| Threshold / categorical | Yes (Gibbs) | Yes | limited | binary/count only | — (Gaussian only) | Yes (`family=`) | Yes (ordinal probit, censored) | Yes (many families + censoring) | — (Gaussian only) |
| Large sparse data | Yes (BLUPF90) | Yes (a core strength) | moderate (dense solver) | Yes (a core strength) | Yes (sparse AI) | single trials | moderate (dense marker matrix) | moderate (MCMC cost) | Yes (big genotype files) |
| Self-contained install | external binaries | Yes | Yes (pure R) | Yes | Yes | Yes (pure R) | Yes (compiled C/Fortran) | Yes | Yes |

lme4 is the general-purpose R workhorse, not a quantitative-genetics tool. Its
sparse-Cholesky engine is fast and its `lmer`/`glmer` interface is the one most
R users already know, but random effects are restricted to i.i.d. grouping
factors: there is no way to supply a covariance matrix for a random term, so an
animal model is out of reach without an extension such as `pedigreemm` (pedigree
annotations) or `lme4qtl` (arbitrary covariance structures). No spatial
correlation, no multi-trait, no genomic features. Reach for it when the design
is a clean nested/crossed hierarchy and you need speed and familiarity.

gremlin is the pure-R REML animal model. It maximises the restricted likelihood
with the average information algorithm and specifies random-effect covariance
through *generalised inverse* matrices, so an **A**⁻¹ from `nadiv` or
`pedigreeTools` drops straight in — as would a **G**⁻¹ or **H**⁻¹ — and sparse
matrix techniques keep it efficient on large pedigrees. The constraints are
real: univariate and Gaussian only, with no genomic, spatial or threshold
machinery. Its appeal is fitting the classical animal model with no external
binaries and no commercial licence.

SpATS solves one problem and solves it well: spatial trend in a field trial,
modelled as a two-dimensional P-spline surface with variance components by REML
via an extension of the SAP (separation of anisotropic penalties) algorithm.
Genotype can be fixed or random (`genotype.as.random`), `family=` gives GLMM
responses through a Fisher-scoring outer loop, and spatial trend and genetic
effects are fitted simultaneously rather than in two stages. It is the de facto
standard for spatial correction in plant breeding, and is often used as a
pre-processing step whose adjusted means feed a downstream genetic analysis.
No relationship matrices, no multi-trait.

MCMCglmm is the closest Bayesian analogue to breedR. It fits multi-response
GLMMs by MCMC, with native animal-model support via a `pedigree=` argument and
arbitrary structured random effects via `ginverse=` — which is also the route to
genomic or single-step relationship matrices, and to spatial (e.g. CAR)
structures. Its distribution list is unusually broad (gaussian, poisson,
categorical, ordinal, threshold, multinomial, zero-inflated and censored
variants), and traits in a multi-response model may follow different
distributions. The costs are the usual MCMC ones: prior specification matters,
chains need diagnosing, and runtimes grow quickly with data size.

BGLR is the odd one out: rather than a general mixed-model fitter, it is a
whole-genome regression package built around a Gibbs sampler. Models are
specified as a list of linear predictor terms (`ETA`), each with its own prior —
`FIXED`, `BRR` (ridge), `BayesA`, `BayesB`, `BayesC`, `BL` (Bayesian LASSO) or
`RKHS`. Pedigree and genomic relationship matrices enter as `RKHS` kernels
(`K = A`, `G` or `H`), which is also how you would fit a spatial or single-step
model — there is no formula interface for pedigrees, AR structures or splines.
The variable-selection priors give posterior inclusion probabilities per marker,
so association work is done through `BFDR()` rather than a dedicated GWAS step.
Everything is posterior samples, not REML point estimates with standard errors.

hibayes overlaps breedR's genomic territory more directly than BGLR does. It
fits Bayesian regressions from three kinds of input: individual-level genotypes
(`ibrm()`), GWAS summary statistics plus an LD matrix (`sbrm()`), and — the
reason it earns a column — single-step models combining genotyped and
ungenotyped animals, with the pedigree passed to `ssbrm()` as an id/sire/dam
matrix. Priors span BayesRR (ridge/GBLUP), A, B, Bpi, C, Cpi, L, R and BSLMM,
and fixed and random effects use ordinary formula syntax with `(1|factor)`.
Association output comes back as posterior inclusion probabilities and window
posterior probability of association, so prediction and GWAS fall out of the
same fit. Single-trait and Gaussian only.

**In short:** ASReml-R is the mature commercial standard, especially for very
large and spatially structured models; sommer is a flexible pure-R package
strong in multivariate genomic prediction; lme4 is the fast, familiar choice for
plain hierarchical designs but has no genetic machinery; gremlin is the
no-dependencies REML animal model, at the price of being univariate and
Gaussian; SpATS is the specialist for spatial correction of field trials;
MCMCglmm is the Bayesian generalist, excellent for multi-response and
non-Gaussian animal models; BGLR is the reference for Bayesian whole-genome
regression and shrinkage/variable-selection priors, but leaves pedigree, spatial
and single-step structures to the user as kernels; hibayes is the closest
genomic competitor, with single-step and GWAS from one Bayesian fit; breedR
combines REML and Bayesian inference with dedicated single-step GBLUP and GWAS,
and competition models, on the scalable BLUPF90 backend at no cost, in exchange
for depending on external binaries.

### See also

Narrower or adjacent tools that did not warrant a column:

- **Mixed models with relationship matrices** — [lme4breeding](https://cran.r-project.org/package=lme4breeding)
  (lme4 plus `relmat`/`addmat` arguments and factor-analytic structures),
  [pedigreemm](https://cran.r-project.org/package=pedigreemm) (lme4 with
  pedigrees), [lme4qtl](https://github.com/variani/lme4qtl) (lme4 with arbitrary
  covariance matrices), [regress](https://cran.r-project.org/package=regress),
  [EMMREML](https://cran.r-project.org/package=EMMREML),
  [rrBLUP](https://cran.r-project.org/package=rrBLUP) (`mixed.solve()`,
  `kin.blup()`, `GWAS()`).
- **Genomic prediction** — [qgg](https://cran.r-project.org/package=qgg)
  (GBLUP, Bayesian regression, marker-set tests),
  [bWGR](https://cran.r-project.org/package=bWGR) (fast multivariate whole-genome
  regression), [SFSI](https://cran.r-project.org/package=SFSI) (sparse selection
  indices), [MegaLMM](https://github.com/deruncie/MegaLMM) (large-scale
  multi-trait), [brms](https://cran.r-project.org/package=brms) (Stan; animal
  models via `gr(id, cov = A)`).
- **GWAS** — [gaston](https://cran.r-project.org/package=gaston),
  [rMVP](https://cran.r-project.org/package=rMVP),
  [GAPIT](https://github.com/jiabowang/GAPIT),
  [statgenGWAS](https://cran.r-project.org/package=statgenGWAS).
- **Spatial and field-trial analysis** — [LMMsolver](https://cran.r-project.org/package=LMMsolver)
  (sparse solver that also implements the SpATS model),
  [statgenSTA](https://cran.r-project.org/package=statgenSTA),
  [nlme](https://cran.r-project.org/package=nlme) (`corStruct` correlated
  errors), [glmmTMB](https://cran.r-project.org/package=glmmTMB) (ar1, ou, exp,
  gau, mat covariance structures).
- **Relationship and pedigree matrices** (inputs, not fitters) —
  [AGHmatrix](https://cran.r-project.org/package=AGHmatrix) (**A**, **D**, **G**,
  **H**, including autopolyploids), [nadiv](https://cran.r-project.org/package=nadiv)
  (dominance and epistatic matrices, inverses),
  [pedigreeTools](https://cran.r-project.org/package=pedigreeTools),
  [ASRgenomics](https://cran.r-project.org/package=ASRgenomics) (**G**-matrix QC
  and tuning).
- **Multi-environment trials** — see
  [ENVIROTYPING.md](ENVIROTYPING.md) for reaction-norm and G×E models driven by
  environmental covariates, plus [metan](https://cran.r-project.org/package=metan)
  and [statgenGxE](https://cran.r-project.org/package=statgenGxE) for AMMI/GGE
  stability analysis.

Sources: [BLUPF90 wiki](https://nce.ads.uga.edu/wiki/), [ASReml-R (VSNi)](https://vsni.co.uk/software/asreml-r/), [sommer (CRAN)](https://cran.r-project.org/package=sommer), [lme4 (CRAN)](https://cran.r-project.org/package=lme4), [Bates *et al.* (2015), *J. Stat. Softw.* 67(1)](https://doi.org/10.18637/jss.v067.i01), [MCMCglmm (CRAN)](https://cran.r-project.org/package=MCMCglmm), [Hadfield (2010), *J. Stat. Softw.* 33(2)](https://doi.org/10.18637/jss.v033.i02), [BGLR (CRAN)](https://cran.r-project.org/package=BGLR), [Pérez & de los Campos (2014), *Genetics* 198:483–495](https://doi.org/10.1534/genetics.114.164442), [Pérez-Rodríguez & de los Campos (2022), *Genetics* 222:iyac112](https://doi.org/10.1093/genetics/iyac112), [gremlin (CRAN)](https://cran.r-project.org/package=gremlin), [SpATS (CRAN)](https://cran.r-project.org/package=SpATS), [Rodríguez-Álvarez *et al.* (2018), *Spat. Stat.* 23:52–71](https://doi.org/10.1016/j.spasta.2017.10.003), [hibayes (CRAN)](https://cran.r-project.org/package=hibayes).

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
