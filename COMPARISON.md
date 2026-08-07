# Comparison with other mixed-model tools

Companion page to the [breedR README](README.md). For packages that bring
environmental covariates into the model, see [ENVIROTYPING.md](ENVIROTYPING.md).

An at-a-glance comparison with eight widely used alternatives for
quantitative-genetics mixed models. They are split by how variance components
are estimated, since that shapes everything else about how a package is used;
breedR appears in both tables because it offers both. Capabilities of all nine
evolve — verify the details against the current documentation for your use case.

## Likelihood-based (REML)

| | **breedR** | **ASReml-R** | **sommer** | **lme4** | **gremlin** | **SpATS** |
|---|---|---|---|---|---|---|
| Engine | BLUPF90 (Fortran) | proprietary (VSNi) | R / C++ (Armadillo) | R / C++ (Eigen) | R / C++ (Matrix) | R (SAP algorithm) |
| License | GPL-3 (free) | commercial | GPL (free) | GPL (≥ 2) | GPL-3 (free) | GPL-2/3 (free) |
| REML algorithm | AI / EM | AI | Newton-Raphson / AI | profiled REML criterion | AI | SAP |
| Pedigree / **A** matrix | Yes (built-in) | Yes | Yes | — (needs `pedigreemm`) | Yes (generalised inverse) | — |
| Spatial (AR, splines) | Yes | Yes | Yes | — | via user-supplied inverse | Yes (2D P-splines) |
| Multi-trait | Yes | Yes | Yes | — | — (univariate only) | — |
| Single-step GBLUP | Yes (dedicated pipeline) | via user-supplied H | via user-supplied H | — | via user-supplied H⁻¹ | — |
| GWAS / SNP effects | Yes (PostGSF90) | via SNP models | Yes | — | — | — |
| Competition / IGE | Yes (built-in) | via custom structures | — | — | — | — |
| Non-Gaussian / threshold | Yes (via Gibbs) | Yes | limited | binary/count only | — (Gaussian only) | Yes (`family=`) |
| Large sparse data | Yes (BLUPF90) | Yes (a core strength) | moderate (dense solver) | Yes (a core strength) | Yes (sparse AI) | single trials |
| Self-contained install | external binaries | Yes | Yes (pure R) | Yes | Yes | Yes (pure R) |

## Bayesian (MCMC)

| | **breedR** | **BGLR** | **MCMCglmm** | **hibayes** |
|---|---|---|---|---|
| Engine | GIBBSF90+ (Fortran) | R / C (Gibbs sampler) | R / C (MCMC) | R / C++ (MCMC) |
| License | GPL-3 (free) | GPL-3 (free) | GPL (≥ 2) | GPL-3 (free) |
| REML also available | Yes | — | — | — |
| Model / prior classes | variance-component Gibbs | `FIXED`, `BRR`, `BayesA/B/C`, `BL`, `RKHS` | Gaussian random effects, inverse-Wishart priors | `BayesRR/A/B/Bpi/C/Cpi/L/R`, `BSLMM` |
| Pedigree / **A** matrix | Yes (built-in) | via RKHS kernel | Yes (`pedigree=`) | Yes (`ssbrm()`) |
| Spatial (AR, splines) | Yes | via user-supplied kernel | via `ginverse=` | — |
| Multi-trait | Yes | Yes (`Multitrait()`) | Yes (a core strength) | — |
| Single-step GBLUP | Yes (dedicated pipeline) | via user-supplied H kernel | via `ginverse=` (H⁻¹) | Yes (`ssbrm()`) |
| GWAS / SNP effects | Yes (PostGSF90) | Yes (BayesB/C + `BFDR()`) | — | Yes (PIP / WPPA) |
| Summary-statistics input | — | — | — | Yes (`sbrm()` + LD matrix) |
| Non-Gaussian / threshold | Yes | Yes (ordinal probit, censored) | Yes (many families + censoring) | — (Gaussian only) |
| Large data | Yes (BLUPF90) | moderate (dense marker matrix) | moderate (MCMC cost) | Yes (big genotype files) |
| Self-contained install | external binaries | Yes (compiled C/Fortran) | Yes | Yes |

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

## See also

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
