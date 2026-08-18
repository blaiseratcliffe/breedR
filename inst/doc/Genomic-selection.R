## ----setup, include = FALSE---------------------------------------------------
# The code in this vignette drives the BLUPF90 backend, which is downloaded at
# install time and is not available when the vignette is built (e.g. on CRAN).
# Chunks are therefore shown but not evaluated; illustrative output is described
# in the text. Run the code in an R session with the binaries installed.
knitr::opts_chunk$set(eval = FALSE, collapse = TRUE, comment = "#>")

## ----install-binaries---------------------------------------------------------
# library(breedR)
# 
# install_genomic_programs()   # preGSf90, postGSf90, predf90, gibbsf90+, ...
# install_renumf90()           # renumf90

## ----data---------------------------------------------------------------------
# data(globulus)
# str(globulus[, 1:7])

## ----genotypes----------------------------------------------------------------
# # A synthetic 200-marker panel for the first 60 animals in the trial.
# set.seed(1)
# gen_ids <- globulus$self[1:60]
# G <- matrix(sample(0:2, length(gen_ids) * 200, replace = TRUE,
#                    prob = c(0.25, 0.5, 0.25)),
#             nrow = length(gen_ids))
# 
# write_snp_file(G, ids = gen_ids, file = "genotypes.txt")
# 
# # Fractional (imputed) dosages are also supported:
# # write_snp_file(dosages, ids = gen_ids, file = "genotypes.txt", fractional = TRUE)

## ----ssgblup------------------------------------------------------------------
# res <- remlf90(
#   fixed   = phe_X ~ gg,
#   genetic = list(model    = 'add_animal',
#                  pedigree = globulus[, 1:3],
#                  id       = 'self'),
#   data    = globulus,
#   genomic = list(snp_file      = "genotypes.txt",
#                  # map_file    = "map.txt",   # optional, for GWAS positions
#                  minfreq       = 0.05,        # MAF filter
#                  save_ginverse = TRUE)        # needed for the GWAS below
# )
# 
# summary(res)

## ----ssgblup-qc---------------------------------------------------------------
# res$genomic$n_snp_total    # markers read
# res$genomic$n_snp_passed   # markers passing QC
# res$genomic$freq           # allele frequencies after cleaning

## ----gwas---------------------------------------------------------------------
# gwas <- postgsf90(res,
#                   windows_variance = 20,     # variance explained per 20-SNP window
#                   manhattan_plot   = TRUE)
# 
# head(gwas$snp_sol)    # per-SNP solutions, weights and variances
# head(gwas$manhattan)  # one row per marker: effect, position (if a map was given)
# head(gwas$windows)    # variance explained by genomic windows

## ----prediction---------------------------------------------------------------
# # genotypes for the animals to predict, in the same BLUPF90 format
# predictions <- predf90(snp_file   = "new_genotypes.txt",
#                        use_mu_hat = TRUE)   # add the base so DGV are on the GEBV scale
# 
# head(predictions)   # data.frame(id, call_rate, dgv)

## ----gibbs--------------------------------------------------------------------
# res_b <- gibbsf90(phe_X ~ gg,
#                   genetic   = list(model = 'add_animal',
#                                    pedigree = globulus[, 1:3], id = 'self'),
#                   data      = globulus,
#                   n_samples = 50000, burnin = 5000, thin = 10)
# 
# colMeans(res_b$samples)              # posterior means of the variance components
# diag <- postgibbsf90(res_b)          # convergence diagnostics (Geweke, ESS, HPD)

