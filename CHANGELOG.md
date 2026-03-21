# Changelog

All notable changes from the original [famuvie/breedR](https://github.com/famuvie/breedR) repository.

## [Unreleased] - 2026-03-20

### BLUPF90+ Migration

The package now uses **BLUPF90+** (the unified program) instead of the deprecated standalone `remlf90` and `airemlf90` binaries. All binaries are downloaded from the official UGA distribution at `nce.ads.uga.edu`.

- **`R/binaries.R`**: Complete rewrite. Switched from breedR GitHub Pages repo (`famuvie.github.io/breedR/bin`) to UGA (`nce.ads.uga.edu/html/projects/programs`). `progsf90_files()` now returns `blupf90+` instead of `c("remlf90", "airemlf90")`. `install_progsf90()` downloads raw executables directly (no archive extraction). Removed `breedR_online()`, `retrieve_bin()`, `decompress()`, `my_unzip()`, `getrootdir()`, `getdir()` — all replaced by `retrieve_bin_direct()`.
- **`R/remlf90-class.R`**: Removed `switch(method, ai=..., em=...)` binary selection. Now uses a single `blupf90+` binary with `OPTION method VCE` and `OPTION EM-REML` for method control.
- **`R/utils.R`**: Simplified `breedR.bin.builtin()` — flat `bin/` directory, no platform-specific subdirectories.
- **`src/install.libs.R`**: Updated for direct download from UGA.
- Supports **Windows**, **Linux**, and **macOS** (via `Mac_OSX/64bit/` path).

### Genomic Selection (ssGBLUP)

Full single-step genomic BLUP pipeline via PREGSF90 integration.

- **New parameter `genomic` in `remlf90()`**: Accepts a list with `snp_file`, `map_file`, QC parameters (`minfreq`, `callrate`, `callrateAnim`, `verify_parentage`), G matrix construction options (`whichG`, `tunedG`, `AlphaBeta`), and save options.
- **Pipeline**: When `genomic` is specified, PREGSF90 runs first (QC + G matrix construction), produces `GimA22i`, then BLUPF90+ runs with `OPTION readGimA22i` for ssGBLUP.
- **QC reports** attached to result as `res$genomic` with allele frequencies, excluded SNPs/animals, Mendelian conflict reports, and summary statistics.

### Genome-Wide Association (ssGWAS)

PostGSF90 wrapper for extracting SNP effects and running GWAS.

- **New function `postgsf90()`**: Takes a fitted genomic `remlf90` object, runs PostGSF90, returns SNP solutions, Manhattan plot data, variance explained by genomic windows, and p-values.
- Supports `windows_variance`, `windows_variance_mbp`, `which_weight`, `manhattan_plot`, `snp_p_value`, and `postgs_trt_eff` options.
- Parses all PostGSF90 output files: `snp_sol`, `chrsnp`, `chrsnp_pval`, `chrsnpvar`, `windows_variance`, `windows_segment`, `snp_pred`.

### Genomic Prediction for New Animals

PREDF90 wrapper for predicting DGV in young/new animals.

- **New function `predf90()`**: Predicts direct genomic values using SNP effects from PostGSF90. Supports `--use_mu_hat`, `--acc`, `--acc_type`, `--use_diagG_acc` options. Returns a data.frame with ID, call rate, DGV, and optionally reliability.

### RENUMF90 Data Preparation

Full wrapper for the RENUMF90 data preparation program.

- **New function `renumf90()`**: Runs RENUMF90 for data renumbering, pedigree validation, unknown parent groups, and inbreeding computation. Accepts either a pre-written parameter file or builds one from R arguments. Parses all output files (`renf90.par`, `renf90.dat`, `renaddxx.ped`, `renf90.tables`, `renf90.inb`).
- **New function `build_renumf90_parfile()`**: Generates RENUMF90 keyword-format parameter files from R arguments. Supports all RENUMF90 keywords including EFFECT, RANDOM, OPTIONAL (pe/mat/mpe), FILE, FILE_POS, SNP_FILE, PED_DEPTH, UPG_TYPE, INBREEDING, and (CO)VARIANCES.
- **New function `remlf90_from_renum()`**: Bridge function that takes RENUMF90 output and runs BLUPF90+ directly, returning solutions and variance components.
- **New functions `install_renumf90()`, `check_renumf90()`**: Binary management for the renumf90 program.

### Bayesian Inference (Gibbs Sampling)

GIBBSF90+ and POSTGIBBSF90 wrappers for Bayesian variance component estimation.

- **New function `gibbsf90()`**: Fits linear mixed models via Gibbs sampling. Same model specification as `remlf90()` (formula, genetic, spatial, genomic). Supports threshold/categorical traits via `cat` parameter (binary survival, ordinal disease scores). Returns posterior samples for variance components, solutions with posterior means and SDs, and deviance for DIC.
- **New function `postgibbsf90()`**: Post-processes Gibbs samples for convergence diagnostics. Returns posterior means, SDs, 95% HPD intervals, effective sample sizes, Geweke diagnostics, autocorrelations at lag 1/10/50, and log marginal density for Bayes factor.
- **New function `gibbsf90_from_renum()`**: Runs GIBBSF90+ on RENUMF90-prepared data. Supports threshold traits and all Gibbs options.
- **New function `build_gibbs_options()`**: Generates GIBBSF90+ OPTION strings (cat, thresholds, prior, seed, fixed_var, solution, cont, save_halfway_samples, hetres_int, residual).
- **Shared pipeline**: Extracted `run_pregsf90_pipeline()` from `remlf90()` into `R/genomic.R` so both `remlf90()` and `gibbsf90()` share the genomic preprocessing code without duplication.

### Genotype Quality Control

- **New function `qcf90()`**: Standalone genotype and pedigree QC using QCF90. Works with non-renumbered (raw alphanumeric) data before RENUMF90. Supports call rate, MAF, monomorphic, HWE, Mendelian conflicts, duplicate detection, chromosome filtering, and clean file output. 57 unit tests.

### Convenience Data Preparation

- **New function `renumf90_from_data()`**: Takes R data.frames and column names instead of file positions. Auto-detects factor vs numeric types, writes temp files, maps column positions, and runs RENUMF90. Validated with globulus dataset. 10 bugs found and fixed during audit.
- **New function `expand_variance()`**: Internal helper for variance matrix expansion (scalar, vector, or matrix to n_traits x n_traits).

### SNP Data Formatting

Functions to read and write genotype files in BLUPF90 format.

- **`write_snp_file()`**: Complete rewrite. Now exported. Supports both integer genotypes (0/1/2/5) and fractional genotypes from imputation (X.XX format). Handles ID padding for fixed-width column alignment. Replaces `NA` with missing value code.
- **New function `read_snp_file()`**: Reads BLUPF90-format genotype files back into R as `list(ids, geno)`. Supports both integer and fractional formats.
- **`write_xref_file()`**: Writes cross-reference ID files for PREGSF90.

### Parallel AR Grid Search

- **New parameter `parallel` in `remlf90()`**: When `parallel = TRUE` or `parallel = <n_cores>`, the AR rho grid search runs BLUPF90+ binaries in parallel using `parallel::parLapply` with socket clusters. Works on all platforms including Windows. Workers only call `system2()` — no R package loading needed.
- Cluster workers are guaranteed to be stopped on failure via `on.exit()` and `tryCatch`. Falls back to sequential with a warning if parallel execution fails.

### Performance Optimizations

- **Vectorized grid indexing** (`R/spatial.R`): Replaced `matrix2vec()` which allocated a full `nx*ny` matrix and used `apply()` per row with direct column-major arithmetic. ~100x faster for 100x100 grids.
- **`rowSums()` instead of `apply(..., 1, sum)`** (`R/simulation.R`): Two call sites in `breedR.sample.phenotype()` replaced with C-level `rowSums()`.

### Error Handling Improvements

- **Binary crash reporting** (`R/remlf90-class.R`): Replaced opaque `stopifnot(is.null(attr(reml.out, 'status')))` with an error message that includes the exit code and last 20 lines of BLUPF90+ output.
- **`renderpf90.default` stops instead of warning** (`R/renderpf90.R`): Unsupported effect types now produce an immediate error listing all supported types, preventing confusing downstream failures.
- **AR grid search resilience** (`R/remlf90-class.R`): Individual rho combinations are wrapped in `tryCatch`. A single failure no longer aborts the entire grid search. Failed combinations get `NA` log-likelihood and are skipped.

### Code Quality

- **Named constants** (`R/breedREnv.R`): Extracted magic numbers to `MAX_COMPETITORS` (8), `MAX_REML_ITERATIONS` (5000), `SPLINE_BOUNDARY_KNOTS` (3), `REGULAR_GRID_QUANTILE_RANGE` (c(.1, .6)).
- **Dead conditional removed** (`R/progsf90.R`): `ifelse(TRUE, 'SE', 'S.D.')` replaced with `'SE'`.
- **Dead code removed**: ~50 lines of commented-out code across `progsf90.R`, `remlf90-class.R`, `splines.R`. Duplicate `match.arg()` call removed from `splines.R`. Unused `build.splines1d.sparse()` function removed. `breedR.result()` test helper moved from production code to `tests/testthat/helper-testdata.R`.
- **`renderpf90.blocks` delegates to `renderpf90.diagonal`** (`R/renderpf90.R`): Eliminated 12 lines of duplicated sparse matrix index extraction code.
- **`renderpf90` contract documented** (`R/AllGeneric.R`): Three-tier return structure (fixed, random, effect_group) fully documented in roxygen with field names, types, and when each tier applies.
- **Pedigree docs fixed** (`R/pedigree.R`): Swapped `@param self` and `@param sire` descriptions which were reversed.

### CI/CD

- **GitHub Actions** (`.github/workflows/R-CMD-check.yaml`): New workflow replacing Travis CI and AppVeyor. 6-cell matrix (ubuntu/windows/macos x R release/oldrel-1). Downloads BLUPF90+ binaries from UGA. INLA install wrapped in `tryCatch` (optional). Integration tests run on ubuntu after R-CMD-check passes.
- **Test coverage** (`.github/workflows/test-coverage.yaml`): Coverage reporting via `covr` and Codecov upload.

### New Tests

- **`tests/testthat/test-refactoring.R`**: 16 tests covering `renderpf90.default` error behavior, `renderpf90.blocks` delegation, vectorized `build_grid` indexing (full grid, sparse grid, column-major ordering), and named constants.

### Prediction Validation

ValidationF90 wrapper for the LR validation method (Legarra & Reverter, 2018).

- **New function `validationf90()`**: Low-level wrapper. Takes parameter file, whole/partial solutions files, and validation IDs. Returns bias, dispersion, and accuracy statistics.
- **New function `validate_prediction()`**: High-level convenience. Automatically fits the model twice (whole and partial data), then runs ValidationF90. Handles partial data creation by setting validation animals' phenotypes to missing.

### Parentage Verification

SeekParentF90 wrapper for SNP-based parentage testing and parent discovery.

- **New function `seekparentf90()`**: Checks parent-offspring Mendelian conflicts, searches for correct parents among candidates, supports multiple SNP chips, year-of-birth filtering, duplicate detection, and trio-based conflict detection.

### Variance Function Helpers

User-friendly functions that generate BLUPF90+ `OPTION se_covar_function` strings for heritability, genetic correlations, and variance proportions with standard errors.

- **New function `h2_formula()`**: Heritability with SE. Handles maternal effects, multi-trait, and additional random effects.
- **New function `rg_formula()`**: Genetic correlation between traits with SE.
- **New function `vp_formula()`**: Total phenotypic variance with SE.
- **New function `var_ratio_formula()`**: Proportion of variance for any random effect.
- **New function `var_functions()`**: Generates all of the above at once for a full model.

### Heterogeneous Residual Variances

- **New function `hetres_options()`**: Generates BLUPF90+ OPTION strings for heterogeneous residual variances modeled as log-linear functions of covariates. Supports both class-based (`hetres_int` for GIBBS3F90) and covariate-based (`hetres_pos`/`hetres_pol` for AIREMLF90) approaches.
- **Known limitation**: `parse_results()` does not yet handle heterogeneous residual output from BLUPF90+ (see issue #1).

### Documentation

- **`inst/doc/pregsf90_postgsf90_options.md`**: Complete reference for all PREGSF90 and PostGSF90 options.
- **`inst/doc/multi_trait_different_effects.md`**: Full workflow guide for multi-trait models with different effects per trait using RENUMF90.
- **`.gitignore`**: Added `inst/bin/`, `.Rhistory`, `.RData`, `.Rproj.user`.

### Files Added

| File | Purpose |
|---|---|
| `R/genomic.R` | Genomic selection: PREGSF90, PostGSF90, PREDF90 wrappers, SNP file I/O, `run_pregsf90_pipeline()` |
| `R/gibbs.R` | GIBBSF90+ and POSTGIBBSF90 wrappers, `gibbsf90()`, `postgibbsf90()`, `gibbsf90_from_renum()` |
| `R/renumf90.R` | RENUMF90 wrapper and `remlf90_from_renum()` bridge |
| `R/validation.R` | ValidationF90 wrapper and `validate_prediction()` convenience |
| `R/seekparent.R` | SeekParentF90 wrapper |
| `R/varfunctions.R` | `h2_formula()`, `rg_formula()`, `var_functions()`, etc. |
| `R/hetres.R` | `hetres_options()` for heterogeneous residual variances |
| `R/qc.R` | QCF90 wrapper for genotype/pedigree quality control |
| `R/renumf90_convenience.R` | `renumf90_from_data()` convenience wrapper |
| `.github/workflows/R-CMD-check.yaml` | GitHub Actions CI |
| `.github/workflows/test-coverage.yaml` | Code coverage workflow |
| `tests/testthat/test-refactoring.R` | Tests for refactored code |
| `tests/testthat/test-variance-functions.R` | 38 tests for variance function helpers |
| `tests/testthat/test-hetres.R` | 13 tests for heterogeneous residual helpers |
| `tests/testthat/test-genomic-helpers.R` | 31 tests for genomic validation/options |
| `tests/testthat/test-snp-io.R` | 20 tests for SNP file read/write roundtrips |
| `tests/testthat/test-renumf90-parfile.R` | 35 tests for RENUMF90 parameter file builder |
| `tests/testthat/test-gibbs-helpers.R` | 28 tests for Gibbs option builder |
| `tests/testthat/test-qcf90.R` | 57 tests for QCF90 argument building and validation |
| `inst/doc/pregsf90_postgsf90_options.md` | PREGSF90/PostGSF90 option reference |
| `inst/doc/multi_trait_different_effects.md` | Multi-trait workflow guide |
| `CHANGELOG.md` | This file |
| `.gitignore` | Git ignore rules |

### Files Modified

| File | Changes |
|---|---|
| `R/binaries.R` | BLUPF90+ migration, UGA source, genomic/renumf90/seekparent/validation binary management |
| `R/remlf90-class.R` | `genomic` and `parallel` params, BLUPF90+ options, PREGSF90 pipeline, parallel grid search, error reporting, hetres docs |
| `R/renderpf90.R` | `default` stops, `blocks` delegates, constants |
| `R/progsf90.R` | Dead code removal, `sd.label` fix, constants |
| `R/spatial.R` | Vectorized indexing, constants |
| `R/simulation.R` | `rowSums()` optimization |
| `R/splines.R` | Dead code removal |
| `R/utils.R` | Simplified `breedR.bin.builtin()`, moved test helper |
| `R/breedREnv.R` | Named constants |
| `R/AllGeneric.R` | `renderpf90` contract documentation |
| `R/pedigree.R` | Swapped `@param` descriptions |
| `src/install.libs.R` | BLUPF90+ migration |
| `tests/testthat/helper-testdata.R` | Added `breedR.result()` test helper |
| `README.md` | Complete rewrite with all new features |
| `DESCRIPTION` | RoxygenNote version bump |
| `NAMESPACE` | New exports |
