# breedR

[![Latest release](https://img.shields.io/github/v/release/blaiseratcliffe/breedR)](https://github.com/blaiseratcliffe/breedR/releases/latest)

### Statistical methods for genetic resources analysts, with genomic selection

breedR is an R package for fitting linear mixed models in breeding and quantitative genetics. It wraps the [BLUPF90](https://nce.ads.uga.edu/wiki/doku.php?id=start) family of programs to estimate variance components by REML, predict breeding values and run genomic evaluations. This fork of the [original breedR](https://github.com/famuvie/breedR) moves it to the current BLUPF90+ program and adds single-step genomic BLUP, GWAS, genomic prediction, Bayesian inference by Gibbs sampling, RENUMF90 data preparation, genotype QC, parentage verification and prediction validation.

> **Note:** This fork is in active development and is written with the assistance of [Claude Code](https://claude.com/claude-code). The genomic paths (ssGBLUP, GWAS, genomic prediction, RENUMF90, parentage verification, and QC) now run end to end and are verified on small test datasets. See [Known limitations](#known-limitations) for two features that need a newer binary or extra setup.
>
> Testing help is welcome. [CONTRIBUTING.md](CONTRIBUTING.md) explains how to compare breedR's estimates with other tools and what to put in an issue.

## Upgrading from an earlier version

Some fixes in 0.13 change results that earlier versions returned without an error. Read "Check these if you have results from an earlier version" in the [v0.13 release notes](https://github.com/blaiseratcliffe/breedR/releases/tag/v0.13), and refit any model that one of its items applies to.

## What breedR fits

| Model | Description |
|---|---|
| Animal model | Additive genetic effects via pedigree |
| Competition model | Direct and competition genetic effects with neighbour structure |
| AR(1)xAR(1) spatial | Autoregressive spatial correlation (row x column) |
| B-splines spatial | Bidimensional penalized splines |
| Blocks | Block or group random effects |
| Generic | User-supplied incidence and covariance or precision matrices |
| Multi-trait | Multiple correlated traits |
| Trait-specific random effects | Random effect groups fitted on some responses of a multi-trait model only (`remlf90(traits = ...)`, [#51](https://github.com/blaiseratcliffe/breedR/issues/51)) |
| Heterogeneous residual variances | Residual variance as a log-linear function of covariates (`hetres_options()`) |
| Threshold/categorical | Binary (survival) and ordinal traits via Gibbs sampling |
| ssGBLUP | Single-step genomic BLUP combining pedigree and genomic data |

## Installation

```r
# The 0.13 release
devtools::install_github('blaiseratcliffe/breedR@v0.13')

# The development version
devtools::install_github('blaiseratcliffe/breedR')
```

Or install from source:

```sh
git clone https://github.com/blaiseratcliffe/breedR.git
R CMD INSTALL breedR
```

breedR needs R >= 3.1.2 and depends on `Matrix`, `sp`, `ggplot2`, `pedigree` and `pedigreemm`, among others. The BLUPF90+ binary is downloaded from [UGA](https://nce.ads.uga.edu/html/projects/programs/) at install time. If that download fails, the package still installs, and `install_progsf90()` fetches the binary later. `install_genomic_programs()` adds preGSf90, postGSf90, predf90, validationf90, predictf90, seekparentf90, gibbsf90+, postgibbsf90 and qcf90, and `install_renumf90()` adds renumf90.

On Linux, downloads from UGA fail at the time of writing with an SSL certificate error, tracked in [#116](https://github.com/blaiseratcliffe/breedR/issues/116). Windows is not affected, and macOS is untested. Most Linux binaries also need the Intel MKL runtime: see [docs/linux.md](docs/linux.md).

## A quick example

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

## Where to go next

- [docs/examples.md](docs/examples.md) has worked examples for multi-trait models, spatial grid search, ssGBLUP, GWAS, Gibbs sampling, RENUMF90, QC, parentage and validation, and lists every wrapped BLUPF90 program.
- The vignettes: `vignette("Overview", package = "breedR")` to start, then `Handling-pedigrees`, `Heritability`, `Heterogeneous-variances`, `Missing-values`, `Additive-Genetic-Models-in-Mixed-Populations`, `General-and-Specific-Combining-Abilities` and `Genomic-selection` (the genomic pipeline from genotypes to prediction).
- [COMPARISON.md](COMPARISON.md) compares breedR with ASReml-R, sommer, lme4 and five other mixed-model packages.
- [ENVIROTYPING.md](ENVIROTYPING.md) covers packages that bring environmental covariates into the model.

## Known limitations

- LR validation with `validate_prediction()` or `validationf90()` needs a newer `validationf90` build. The wrapper is verified, and the whole and partial fits and the partial data it builds are correct, but the bundled `validationf90` v1.01 aborts on valid input with an end-of-file read error inside the program.
- `predf90(acc = TRUE)` cannot compute reliabilities yet: they need an `OPTION snp_var` file from the `postgsf90()` run, which breedR does not write. Direct genomic values (`acc = FALSE`, the default) work.

## Development

```bash
# Run the unit tests (loads the package first)
Rscript -e "devtools::test()"

# Regenerate documentation
Rscript -e "roxygen2::roxygenise()"

# Build and check the package
R CMD build . --no-build-vignettes
R CMD check breedR_*.tar.gz --no-vignettes
```

## Credits

breedR was originally developed by [Facundo Munoz](https://github.com/famuvie) as part of the Trees4Future and ProCoGen projects. The BLUPF90 programs are developed by [Ignacy Misztal's group](https://nce.ads.uga.edu/) at the University of Georgia.

This fork is maintained by [Blaise Ratcliffe](https://github.com/blaiseratcliffe).

## License

GPL-3
