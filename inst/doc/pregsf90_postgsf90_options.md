# PreGSF90 / PostGSF90 Complete Option Reference

Source: https://nce.ads.uga.edu/wiki/doku.php?id=readme.pregsf90
Last modified: 2025/10/03

## Input Files

### Required
- **Parameter file** (renf90.par) with `OPTION SNP_file <filename>`
- **Genotype file**: Field 1 = animal ID (same format as pedigree), Field 2 = genotypes (0/1/2/5=missing, fixed-width). Fractional genotypes supported (two decimal places, no spaces).
- **XrefID file** (`<genotype_file>_XrefID`): renumbered ID + original ID, created by RENUMF90

### Optional
- Allele frequencies file (`OPTION FreqFile`)
- Map file (`OPTION map_file`) — header: SNP_ID, CHR, POS
- Weight file (`OPTION weightedG`)
- Pre-computed G, A22, or their inverses (via save/read options)

### Input Format Options
- `OPTION missingAIPL` — read genotype codes 3,4 as missing (converts to 5)
- `OPTION QMSim` — read genotype codes 3,4 as 1 (heterozygote)
- `OPTION fastread` — use C zlib for faster reading (large datasets)
- `OPTION maxsnp x` — max string length for marker data (default 400,000)

## Output Files (PreGSF90)

| File | Content |
|---|---|
| `GimA22i` | inv(G) - inv(A22), binary format (default) |
| `freqdata.count` | Allele frequencies before QC |
| `freqdata.count.after.clean` | Allele frequencies after QC + exclusion codes |
| `Gen_call_rate` | Animals excluded by call rate |
| `Gen_conflicts` | Mendelian conflict report |

### Exclusion Codes (freqdata.count.after.clean)
1. Call Rate
2. MAF
3. Monomorphic
4. Excluded by request
5. Mendelian error
6. HWE
7. High correlation with other SNP

## G Matrix Construction

### Method
| Option | Syntax | Description | Default |
|---|---|---|---|
| whichG | `OPTION whichG x` | 1: ZZ'/k (VanRaden 2008), 2: ZDZ'/n (Amin 2007), 3: as 2 + UAR (Yang 2010) | 1 |
| whichfreq | `OPTION whichfreq x` | 0: from file, 1: 0.5, 2: from genotypes | 2 |
| whichfreqScale | `OPTION whichfreqScale x` | 0: from file, 1: 0.5, 2: from genotypes | 2 |
| FreqFile | `OPTION FreqFile <file>` | Read allele frequencies from file (format: snp, freq) | freqdata |
| whichScale | `OPTION whichScale x` | 1: 2sum(p(1-p)), 2: tr(ZZ')/n, 3: Gianola correction | 1 |
| weightedG | `OPTION weightedG <file>` | Weighted G: Z*=Z*sqrt(D), G=ZDZ' | - |

## Quality Control

### Default QC (always run unless disabled)
- MAF
- Call rate (SNPs and animals)
- Monomorphic
- Parent-progeny conflicts (SNPs and animals)

### QC Parameters
| Option | Syntax | Default | Description |
|---|---|---|---|
| minfreq | `OPTION minfreq x` | 0.05 | Ignore SNP with MAF < x |
| callrate | `OPTION callrate x` | 0.90 | Ignore SNP with call rate < x |
| callrateAnim | `OPTION callrateAnim x` | 0.90 | Ignore animals with call rate < x |
| monomorphic | `OPTION monomorphic x` | 1 (enabled) | 0 to disable |
| hwe | `OPTION hwe x` | NOT run | Max observed-expected freq diff (default 0.15) |
| high_correlation | `OPTION high_correlation x y` | NOT run | x=freq diff threshold (0.025), y=identity threshold (0.995) |
| verify_parentage | `OPTION verify_parentage x` | 3 | 0: none, 1: detect, 2: search alternate (SeekParentF90), 3: detect+eliminate |
| exclusion_threshold | `OPTION exclusion_threshold x` | 1 | % of SNP for wrong relationship |
| exclusion_threshold_snp | `OPTION exclusion_threshold_snp x` | 10 | % exclusions per locus to exclude SNP |
| number_parent_progeny_evaluations | `OPTION number_parent_progeny_evaluations x` | 100 | Min pairs to exclude SNP |
| no_quality_control | `OPTION no_quality_control` | - | Disable ALL QC |
| outcallrate | `OPTION outcallrate` | - | Print call rate files |
| h2_gene_content | `OPTION h2_gene_content` | - | Check h2 of gene content (Forneris 2015) |

### Parentage/Conflict Reports
| Option | Syntax | Description |
|---|---|---|
| outparent_progeny | `OPTION outparent_progeny` | Full log of all parent-progeny pairs tested |
| out_snp_exclusion_error_rate | `OPTION out_snp_exclusion_error_rate` | SNP Mendelian error rate stats |

### Chromosome/Sample Filtering
| Option | Syntax | Description |
|---|---|---|
| excludeCHR | `OPTION excludeCHR n1 n2 ...` | Exclude chromosomes (needs map_file) |
| includeCHR | `OPTION includeCHR n1 n2 ...` | Include only these chromosomes |
| excludeSample | `OPTION excludeSample n1 n2 ...` | Exclude samples by row number |
| sex_chr | `OPTION sex_chr n` | Chromosomes >= n are not autosomes |

### G Matrix Diagnostics
| Option | Syntax | Default | Description |
|---|---|---|---|
| threshold_duplicate_samples | `OPTION threshold_duplicate_samples x` | 0.9 | Warn if G(i,j)/sqrt(G(i,i)*G(j,j)) > x |
| high_threshold_diagonal_g | `OPTION high_threshold_diagonal_g x` | 1.6 | Flag high G diagonals |
| low_threshold_diagonal_g | `OPTION low_threshold_diagonal_g x` | 0.7 | Flag low G diagonals |
| plotpca | `OPTION plotpca <print/noprint>` | - | PCA for stratification |
| extra_info_pca | `OPTION extra_info_pca <file> col` | - | Color PCA by class variable |

### LD Calculation
| Option | Syntax | Default | Description |
|---|---|---|---|
| calculate_LD | `OPTION calculate_LD` | - | Squared correlation of allele counts |
| LD_by_chr | `OPTION LD_by_chr` | - | LD within chromosome |
| LD_by_pos | `OPTION LD_by_pos x` | 200000 bp | LD within chromosome windows |
| filter_by_LD | `OPTION filter_by_LD x` | 0.8 | Filter SNP with Rsq > threshold |
| thr_output_LD | `OPTION thr_output_LD x` | 0.1 | Print threshold for Rsq |

### Clean Data Output
| Option | Syntax | Description |
|---|---|---|
| saveCleanSNPs | `OPTION saveCleanSNPs` | Save QC'd genotype data (gt_clean, gt_clean_XrefID, gt_SNPs_removed, gt_Animals_removed) |

## A22/G Correlation QC
| Option | Syntax | Default | Description |
|---|---|---|---|
| thrWarnCorAG | `OPTION thrWarnCorAG x` | 0.5 | Warn if cor(A22,G) < x |
| thrStopCorAG | `OPTION thrStopCorAG x` | 0.3 | Stop if cor(A22,G) < x |
| thrCorAG | `OPTION thrCorAG x` | 0.02 | Only calc cor for A22 >= x |

## H Matrix Options

Formula: `tau * inv(alpha*G + beta*A22 + gamma*I + delta*11') - omega * inv(A22)`

Defaults: tau=1, alpha=0.95, beta=0.05, gamma=0, delta=0, omega=1

| Option | Syntax |
|---|---|
| TauOmega | `OPTION TauOmega tau omega` |
| AlphaBeta | `OPTION AlphaBeta alpha beta` |
| GammaDelta | `OPTION GammaDelta gamma delta` |
| tunedG | `OPTION tunedG x` — 0: none, 1: diag=1/offdiag=0, 2: match A22 means (default), 3: mean(G)=mean(A22), 4: Fst (Powell 2010), 9: arbitrary `a+b*G` |

## Diagonal of H (Genomic Inbreeding)

| Option | Syntax | Description |
|---|---|---|
| saveDiagH | `OPTION saveDiagH` | diag(H) with renumbered IDs |
| saveDiagHOrig | `OPTION saveDiagHOrig` | diag(H) with original IDs |
| methodDiagH | `OPTION methodDiagH x` | 1: sparse inversion (default, fast for small/medium), 2: outer product (recommended for large pedigrees) |

## GWAS Options (PostGSF90)

### Manhattan Plots
| Option | Syntax | Description |
|---|---|---|
| Manhattan_plot | `OPTION Manhattan_plot` | Via GNUPLOT |
| Manhattan_plot_R | `OPTION Manhattan_plot_R` | Via R (creates PDF/PNG/TIF) |
| Manhattan_plot_R_format | `OPTION Manhattan_plot_R_format <format>` | pdf (default), png, tif |
| plotsnp | `OPTION plotsnp n` | 1: abs(val), 2: abs(val/sd) (default) |

### SNP Analysis
| Option | Syntax | Description |
|---|---|---|
| SNP_moving_average | `OPTION SNP_moving_average n` | Moving average of n adjacent SNPs |
| windows_variance | `OPTION windows_variance n` | Variance by n adjacent SNPs (moving) |
| windows_variance_mbp | `OPTION windows_variance_mbp n` | Variance by n Mb windows |
| windows_variance_type | `OPTION windows_variance_type n` | 1: moving, 2: exclusive windows |
| which_weight | `OPTION which_weight x` | 1: y^2*(2p(1-p)), 2: y^2, 3: experimental, 4/nonlinearA: C**(abs(y)/sqrt(var)-2) |
| solutions_postGS | `OPTION solutions_postGS x` | Solutions filename (default: solutions) |
| postgs_trt_eff | `OPTION postgs_trt_eff x1 x2` | Compute for specific trait x1, effect x2 |
| snp_p_value | `OPTION snp_p_value` | Compute p-values (expensive) |
| snp_var | `OPTION snp_var` | PEC for SNP (for PREDF90 reliability) |

### PostGSF90 Output Files
| File | Columns |
|---|---|
| snp_sol | trait, effect, SNP, Chr, Pos, solution, weight, [variance], [var_solution] |
| chrsnp | trait, effect, abs(SNP_i)/SD(SNP), SNP, Chr, Pos |
| chrsnp_pval | trait, effect, -log10(p-value), SNP, Chr, Pos |
| chrsnpvar | trait, effect, variance, SNP, Chr, Pos |
| windows_segment | label, window_size, start_SNP, end_SNP, window_id, start_pos, end_pos |
| windows_variance | trait, effect, start_SNP, end_SNP, window_size, start_pos, end_pos, window_id, variance |
| snp_pred | allele frequencies + SNP effects |

## Save/Read Intermediate Files

### Save Options
| Option | Description |
|---|---|
| `OPTION saveAscii` | Save matrices as ASCII (default: binary) |
| `OPTION saveHinv` | Save H^-1 in Hinv.txt |
| `OPTION saveAinv` | Save A^-1 in Ainv.txt |
| `OPTION saveHinvOrig` | H^-1 with original IDs |
| `OPTION saveAinvOrig` | A^-1 with original IDs |
| `OPTION saveDiagGOrig` | Diagonal of G with original IDs |
| `OPTION saveGOrig` | G with original IDs |
| `OPTION saveA22Orig` | A22 with original IDs |
| `OPTION saveGimA22iOrig` | GimA22i with original IDs |
| `OPTION saveGimA22iRen` | GimA22i with renumbered IDs |
| `OPTION savePLINK` | Genotypes in PLINK format |
| `OPTION no_full_binary` | Half-matrix (compatibility) |
| `OPTION saveA22` | Save A22 |
| `OPTION saveA22Inverse` | Save A22^-1 |
| `OPTION saveA22InverseOrig` | A22^-1 with original IDs |
| `OPTION saveG [all]` | Save G (all: save all intermediates) |
| `OPTION saveGInverse` | Save G^-1 |
| `OPTION saveGInverseOrig` | G^-1 with original IDs |
| `OPTION saveGmA22` | Save G - A22 |
| `OPTION saveCleanSNPs` | Save QC'd data |
| `OPTION readOrigId` | Read original IDs from renaddxx.ped |

### Read Options
| Option | Description |
|---|---|
| `OPTION readGimA22i [file]` | Read pre-computed GimA22i |
| `OPTION readG [file]` | Read G |
| `OPTION readGInverse [file]` | Read G^-1 |
| `OPTION readA22 [file]` | Read A22 |
| `OPTION readA22Inverse [file]` | Read A22^-1 |
| `OPTION readGmA22 [file]` | Read G - A22 |

## Misc Options
| Option | Syntax | Description |
|---|---|---|
| num_threads_pregs | `OPTION num_threads_pregs n` | MKL-OpenMP threads for matrices |
| num_threads_iod | `OPTION num_threads_iod n` | MKL-OpenMP threads for PCG |
| graphics | `OPTION graphics s` | GNUPLOT plots (s = display seconds) |
| msg | `OPTION msg x` | 0: minimal, 1: verbose |

## Deprecated
| Option | Replacement |
|---|---|
| `OPTION chrinfo file` | Use `OPTION map_file` instead |
