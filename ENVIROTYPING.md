# Modelling environmental covariates (enviromics and reaction norms)

A companion to the [tool comparison](COMPARISON.md), which covers packages that
fit mixed models over pedigrees, markers and space. This page covers a different
problem: bringing *environmental*
information — weather, soil, elevation, management — into the model as
covariates, so that genotype × environment interaction becomes something you can
predict rather than merely describe.

Capabilities evolve quickly in this area and several of these packages live on
GitHub rather than CRAN. Verify the details against current documentation before
committing to a pipeline.

## The reaction-norm framework

Nearly everything below descends from the model of
[Jarquín *et al.* (2014)](https://doi.org/10.1007/s00122-013-2243-1). Given a
matrix **W** of *p* centred and scaled environmental covariates measured on each
environment, define an environmental relationship matrix

    Ω = W W' / p

by analogy with the genomic relationship matrix **G** = XX'/p. The phenotype is
then modelled as

    y = 1μ + Xβ + u_g + u_w + u_gw + ε

with

| Term | Distribution | Interpretation |
|---|---|---|
| `u_g` | N(0, **G** σ²_g) | main effect of the genotype |
| `u_w` | N(0, **Ω** σ²_w) | main effect of the environment |
| `u_gw` | N(0, (**G** ∘ **Ω**) σ²_gw) | reaction norm — the G×E term |

The Hadamard (element-wise) product **G** ∘ **Ω** is the key idea: it replaces an
unstructured G×E term, which cannot extrapolate, with a structured one that
*can* — because a new environment enters the model through its covariates, not
through its identity. That is what makes prediction into untested environments
possible.

The practical consequence is that the work splits into two halves, and most
packages do only one of them:

1. **Build W** — retrieve weather and soil data, align it to planting dates and
   growth stages, compute agroclimatic indices, and summarise to one row per
   environment. This is "envirotyping".
2. **Fit the model** — construct the kernels and estimate variance components.
   Any package that accepts user-supplied covariance matrices can do this,
   including breedR.

## 1. Envirotyping pipelines (building W)

| Package | Where | Notes |
|---|---|---|
| **EnvRtype** | [GitHub](https://github.com/allogamous/EnvRtype) / [docs](https://allogamous.github.io/EnvRtype/) | The reference implementation. Three modules: environmental sensing (`get_weather()`, `processWTH()`), macro-environmental characterisation (`env_typing()`, `W_matrix()`), and prediction (`env_kernel()`, `get_kernel()`, `kernel_model()`). Derives ecophysiological variables (`param_temperature()`, `param_radiation()`, `param_atmospheric()`) rather than just raw weather. Fits models through BGGE. |
| **envirotypeR** | [GitHub](https://github.com/gcostaneto/envirotypeR) | Costa-Neto's successor to EnvRtype (MIT, not on CRAN), rewritten to address user feedback since 2021 and broadened to "plant and animal enviromics". `get_weather()` (NASA POWER), `get_soil()` (SoilGrids), `get_spatial()`, `WC_Bioclimate()` (WorldClim), `SRTM_elevation()`, `process_synthetic()`. |
| **learnMET** | [GitHub](https://cjubin.github.io/learnMET/) | Bundles covariate construction with ML fitting. Retrieves NASA POWER data or accepts field weather stations, and aggregates daily weather over either naive fixed windows (e.g. non-overlapping 10-day periods) or phenological stages — the aggregation question is usually the hard part, and this is one of the few packages to treat it explicitly. |

## 2. Fitting reaction-norm and G×E kernel models

| Package | Approach | Notes |
|---|---|---|
| **BGLR** | Bayesian, MCMC | The original Jarquín implementation. Each kernel enters `ETA` as an `RKHS` term with `K = G`, `Ω` or `G ∘ Ω`. Slow but completely general. |
| **BGGE** | Bayesian, MCMC | Purpose-built successor for exactly this model class: one function to prepare the genomic G×E kernels, one to fit them. Reported accuracies comparable to BGLR at up to five times lower computational cost. The engine behind EnvRtype. |
| **BMTME** | Bayesian, MCMC | Multi-trait multi-environment models, Montesinos-López *et al.* |
| **sommer** / **ASReml-R** / **lme4breeding** | REML | General mixed-model fitters that accept user-supplied covariance matrices, so **Ω** and **G** ∘ **Ω** can be passed directly. Factor-analytic structures are the classical alternative to covariate-driven G×E here. |
| **MegaLMM** | REML / variational | Large-scale multi-trait and multi-environment mixed models. |
| **SFSI** | penalised | Sparse selection indices and multi-trait/environment sparse genomic prediction; borrows information across correlated traits and environments from subsets of the training set. |
| **bWGR** | Bayesian, fast | Whole-genome regression with multivariate/multi-environment models under a unified Bayesian GBLUP framework. |
| **breedR** | REML / Bayesian | See [below](#where-breedr-fits). |

## 3. Reaction norms without explicit covariates

Worth knowing about, because they answer the same biological question with far
less data plumbing — they regress on an *environmental index* derived from the
phenotypes themselves rather than on measured covariates.

| Package | Notes |
|---|---|
| **FW** | Finlay–Wilkinson regression, both the classical two-step OLS version and a Bayesian single-step version that incorporates pedigree or marker-derived kinship for varieties and a covariance structure for environments. Parsimonious, and often a sensible baseline before reaching for full enviromics. |
| **CERIS-JGRA** | Scripts rather than a packaged release. Searches for the window of the growing season in which an environmental variable best explains the environmental index, which addresses the same aggregation problem `learnMET` tackles from the ML side. |

## 4. Machine-learning approaches

| Package | Notes |
|---|---|
| **learnMET** | Gradient-boosted trees, random forests, stacked ensembles and multi-layer perceptrons over combined genomic and environmental features, with cross-validation schemes designed for MET data (predicting new genotypes, new environments, or both). |
| **MTMEGPS** | Deep learning for multi-trait, multi-environment genomic and phenomic selection. |

The appeal is that trees and networks capture non-linear and threshold responses
to environment that a linear reaction norm cannot. The cost is the usual one:
variance components and heritabilities are no longer available, so these
complement rather than replace the mixed-model route.

## 5. Where the environmental data comes from

None of these fit models; they supply the raw material for **W**.

| Package | Source |
|---|---|
| **nasapower** | NASA POWER global meteorology, solar and climatology, ~0.5° resolution. The default backend for both EnvRtype and envirotypeR. |
| **chirps** | CHIRPS high-resolution precipitation. |
| **climatrends** | Computes established temperature, precipitation and crop-sensitive indices from any of the above. Explicitly designed for heterogeneous, decentralised trials where locations differ in climate *and* in planting date and season length. |
| **geodata** | WorldClim bioclimatic variables, SoilGrids, elevation. |
| **daymetr** | Daymet daily surface weather, North America. |
| **easyclimate** | High-resolution daily climate for Europe. |
| **soilDB** | USDA-NRCS soil survey databases. |
| **apsimx** | Interface to the APSIM crop model, including weather-file retrieval — the route to *simulated* stress covariates (e.g. modelled water deficit) rather than raw weather. |

Crop-model-derived covariates deserve a note: a simulated drought-stress index
over a defined growth stage is often a far better predictor than the rainfall
that produced it, because the model does the biological aggregation for you.

## 6. Classical G×E analysis (descriptive, no covariates)

Different goal — partitioning and visualising an existing G×E signal rather than
predicting into new environments — but frequently the first step in the same
analysis.

| Package | Notes |
|---|---|
| **metan** | Broad multi-environment trial analysis: AMMI, GGE biplots, ecovalence, joint regression, and BLUP-based stability statistics including WAASB, plus cross-validation for AMMI and BLUP models. |
| **statgenGxE** | Biometris MET suite: `gxeAmmi()`, `gxeGGE()`, Finlay–Wilkinson, stability measures, with strong plotting and reporting. Pairs with `statgenSTA` for single-trial spatial analysis. |
| **geneticae** | AMMI and GGE with attention to imputation of unbalanced trials. |

## Where breedR fits

breedR has no envirotyping pipeline and no weather retrieval — those belong to
the packages in sections 1 and 5, and are best run as a separate upstream step.

What breedR does provide is the model-fitting half. The `generic()` effect takes
a user-supplied incidence matrix together with a covariance or precision matrix,
which is exactly the interface the reaction-norm model needs:

```r
## G: genotypes x genotypes genomic relationship matrix
## W: environments x covariates, centred and scaled
Omega <- tcrossprod(W) / ncol(W)             # environmental relationship matrix

## Incidence matrices, one row per record
Z_g  <- model.matrix(~ 0 + g,  dat)          # dat$g  = genotype factor
Z_e  <- model.matrix(~ 0 + e,  dat)          # dat$e  = environment factor
Z_ge <- model.matrix(~ 0 + ge, dat)          # dat$ge = genotype:environment

## Hadamard product, expanded to the genotype x environment combinations
GxW <- G[dat$gid, dat$gid] * Omega[dat$eid, dat$eid]

res <- remlf90(
  fixed   = y ~ 1,
  generic = list(
    gen = list(incidence = Z_g,  covariance = G),
    env = list(incidence = Z_e,  covariance = Omega),
    gxe = list(incidence = Z_ge, covariance = GxW)
  ),
  data = dat
)
```

`summary(res)` then reports σ²_g, σ²_w, σ²_gw and σ²_e, and `ranef(res)$gxe`
holds the reaction-norm deviations. Two practical notes: nudge the diagonals
(`diag(K) <- diag(K) + 1e-4`) if a kernel is singular, and expect the G×E term
to be poorly separated from the residual unless there is replication within
environments — with a single record per genotype × environment the two compete
for the same variance.

Because the backend is BLUPF90, this scales to MET data that would be awkward
for a pure-R sampler, and the same kernels can be fitted by REML or — via
`gibbsf90()` — by MCMC. The trade-off relative to EnvRtype/BGGE is that you
assemble the kernels yourself and get no envirotyping conveniences.

## References

- Jarquín D, Crossa J, Lacaze X, Du Cheyron P, Daucourt J, *et al.* (2014).
  A reaction norm model for genomic selection using high-dimensional genomic and
  environmental data. *Theor Appl Genet* 127:595–607.
  [doi:10.1007/s00122-013-2243-1](https://doi.org/10.1007/s00122-013-2243-1)
- Costa-Neto G, Galli G, Carvalho HF, Crossa J, Fritsche-Neto R (2021).
  EnvRtype: a software to interplay enviromics and quantitative genomics in
  agriculture. *G3* 11:jkab040.
  [doi:10.1093/g3journal/jkab040](https://doi.org/10.1093/g3journal/jkab040)
- Granato I, Cuevas J, Luna-Vázquez F, Crossa J, Montesinos-López O, Burgueño J,
  Fritsche-Neto R (2018). BGGE: a new package for genomic-enabled prediction
  incorporating genotype × environment interaction models. *G3* 8:3039–3047.
  [doi:10.1534/g3.118.200435](https://doi.org/10.1534/g3.118.200435)
- Lian L, de los Campos G (2016). FW: an R package for Finlay–Wilkinson
  regression that incorporates genomic/pedigree information and covariance
  structures between environments. *G3* 6:589–597.
  [doi:10.1534/g3.115.026328](https://doi.org/10.1534/g3.115.026328)
- Westhues CC, Simianer H, Beissinger TM (2022). learnMET: an R package to apply
  machine learning methods for genomic prediction using multi-environment trial
  data. *G3* 12:jkac226.
  [doi:10.1093/g3journal/jkac226](https://doi.org/10.1093/g3journal/jkac226)
- Olivoto T, Lúcio AD (2020). metan: an R package for multi-environment trial
  analysis. *Methods Ecol Evol* 11:783–789.
  [doi:10.1111/2041-210X.13384](https://doi.org/10.1111/2041-210X.13384)
- de Sousa K, van Etten J, Solberg SØ (2023). Climate variability indices for
  ecological and crop models in R: the climatrends package. *J Open Source
  Softw*. [doi:10.21105/joss.04405](https://doi.org/10.21105/joss.04405)
