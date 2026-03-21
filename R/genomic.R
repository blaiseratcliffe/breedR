### Genomic selection support functions ###
### Wraps PREGSF90 for single-step GBLUP (ssGBLUP)


#' Validate genomic specification
#'
#' Checks the genomic argument to remlf90() and sets defaults.
#'
#' @param genomic a list with genomic parameters. See \code{\link{remlf90}}.
#' @return The validated and defaulted genomic list.
check_genomic <- function(genomic) {

  if (!is.list(genomic))
    stop("'genomic' must be a list.", call. = FALSE)

  # snp_file is required
  if (is.null(genomic$snp_file))
    stop("'genomic$snp_file' is required.", call. = FALSE)
  if (!file.exists(genomic$snp_file))
    stop("SNP file not found: ", genomic$snp_file, call. = FALSE)

  # map_file is optional but must exist if provided
  if (!is.null(genomic$map_file) && !file.exists(genomic$map_file))
    stop("Map file not found: ", genomic$map_file, call. = FALSE)

  # G matrix method: 1 (VanRaden), 2 (Amin), 3 (Yang)
  if (is.null(genomic$whichG)) genomic$whichG <- 1L
  if (!genomic$whichG %in% 1:3)
    stop("'genomic$whichG' must be 1, 2, or 3.", call. = FALSE)

  # G matrix tuning: 0-4 or 9
  if (is.null(genomic$tunedG)) genomic$tunedG <- 2L
  if (!genomic$tunedG %in% c(0:4, 9))
    stop("'genomic$tunedG' must be 0, 1, 2, 3, 4, or 9.", call. = FALSE)

  # Blending weights
  if (is.null(genomic$AlphaBeta)) genomic$AlphaBeta <- c(0.95, 0.05)
  if (!is.numeric(genomic$AlphaBeta) || length(genomic$AlphaBeta) != 2)
    stop("'genomic$AlphaBeta' must be a numeric vector of length 2.", call. = FALSE)

  # QC thresholds
  if (is.null(genomic$minfreq)) genomic$minfreq <- 0.05
  if (is.null(genomic$callrate)) genomic$callrate <- 0.90
  if (is.null(genomic$callrateAnim)) genomic$callrateAnim <- 0.90

  for (param in c('minfreq', 'callrate', 'callrateAnim')) {
    val <- genomic[[param]]
    if (!is.numeric(val) || length(val) != 1 || val < 0 || val > 1)
      stop("'genomic$", param, "' must be a single number between 0 and 1.",
           call. = FALSE)
  }

  # Parentage verification: 0 (none), 1 (detect), 3 (detect+eliminate)
  if (is.null(genomic$verify_parentage)) genomic$verify_parentage <- 3L
  if (!genomic$verify_parentage %in% 0:3)
    stop("'genomic$verify_parentage' must be 0, 1, 2, or 3.", call. = FALSE)

  # Save options
  if (is.null(genomic$saveG)) genomic$saveG <- FALSE
  if (is.null(genomic$saveA22)) genomic$saveA22 <- FALSE

  # Extra raw options (character vector)
  if (!is.null(genomic$extra_options) && !is.character(genomic$extra_options))
    stop("'genomic$extra_options' must be a character vector.", call. = FALSE)

  return(genomic)
}


#' Build PROGSF90 option lines for genomic analysis
#'
#' Converts a validated genomic specification into a character vector of
#' OPTION lines for the PROGSF90 parameter file.
#'
#' @param genomic a validated genomic list (output of check_genomic).
#' @return Character vector of option strings (without the "OPTION" prefix,
#'   which is added by write.progsf90).
build_genomic_options <- function(genomic) {

  opts <- character(0)

  # Required: SNP file (use basename since files are copied to tmpdir)
  opts <- c(opts, paste("SNP_file", basename(genomic$snp_file)))

  # Map file
  if (!is.null(genomic$map_file))
    opts <- c(opts, paste("map_file", basename(genomic$map_file)))

  # G matrix construction
  opts <- c(opts, paste("whichG", genomic$whichG))
  opts <- c(opts, paste("tunedG", genomic$tunedG))
  opts <- c(opts, paste("AlphaBeta",
                         genomic$AlphaBeta[1], genomic$AlphaBeta[2]))

  # QC thresholds
  opts <- c(opts, paste("minfreq", genomic$minfreq))
  opts <- c(opts, paste("callrate", genomic$callrate))
  opts <- c(opts, paste("callrateAnim", genomic$callrateAnim))
  opts <- c(opts, paste("verify_parentage", genomic$verify_parentage))

  # Save options
  if (isTRUE(genomic$saveG)) opts <- c(opts, "saveG")
  if (isTRUE(genomic$saveA22)) opts <- c(opts, "saveA22")

  # Save ASCII for inspection
  if (isTRUE(genomic$saveAscii)) opts <- c(opts, "saveAscii")

  # Extra raw options passed through directly
  if (!is.null(genomic$extra_options))
    opts <- c(opts, genomic$extra_options)

  return(opts)
}


#' Run the full PREGSF90 genomic preprocessing pipeline
#'
#' Handles file copying, parameter file writing, PREGSF90 execution, and
#' GimA22i verification. Used internally by \code{remlf90()} and
#' \code{gibbsf90()} when \code{genomic} is specified.
#'
#' @param genomic validated genomic specification (from check_genomic).
#' @param genomic_opts character vector of genomic OPTION lines.
#' @param pf90 progsf90 object.
#' @param tmpdir working directory.
#' @param breedR.bin binary directory.
#' @return List with updated \code{pf90} and \code{pregs_out}.
run_pregsf90_pipeline <- function(genomic, genomic_opts, pf90,
                                   tmpdir, breedR.bin) {
  # Write parameter file for PREGSF90 with genomic options only
  pf90_pregs <- pf90
  pf90_pregs$parameter$options <- genomic_opts
  write.progsf90(pf90_pregs, dir = tmpdir)

  # Copy input files to working directory (skip if already there)
  copy_if_needed <- function(src) {
    dest <- file.path(tmpdir, basename(src))
    if (normalizePath(src, mustWork = FALSE) !=
        normalizePath(dest, mustWork = FALSE))
      file.copy(src, dest, overwrite = TRUE)
  }

  copy_if_needed(genomic$snp_file)
  if (!is.null(genomic$map_file))
    copy_if_needed(genomic$map_file)
  xref_file <- paste0(genomic$snp_file, "_XrefID")
  if (file.exists(xref_file))
    copy_if_needed(xref_file)

  # Check binaries
  if (!check_genomic_programs(breedR.bin, quiet = TRUE))
    stop("Genomic program binaries (preGSf90) not installed. ",
         "See ?install_genomic_programs", call. = FALSE)

  # Run PREGSF90
  pregs_out <- run_pregsf90(tmpdir, breedR.bin)
  writeLines(pregs_out, file.path(tmpdir, "pregsf90.out"))

  # Verify GimA22i
  if (!file.exists(file.path(tmpdir, "GimA22i")))
    stop("preGSf90 did not produce GimA22i. Check QC reports:\n",
         file.path(tmpdir, "pregsf90.out"), call. = FALSE)

  # Update pf90 options for downstream program
  pf90$parameter$options <- c(pf90$parameter$options,
                               paste("SNP_file", basename(genomic$snp_file)),
                               "readGimA22i")

  return(list(pf90 = pf90, pregs_out = pregs_out))
}


#' Run PREGSF90 binary
#'
#' Executes the preGSf90 program in the specified directory.
#' The directory must already contain the parameter file, data, pedigree,
#' genotype file, and XrefID file.
#'
#' @param dir working directory containing all input files.
#' @param bin_path directory containing the preGSf90 binary.
#' @return Character vector of program stdout.
run_pregsf90 <- function(dir, bin_path) {

  pregs_name <- genomic_program_files(breedR.os.type())[1]
  pregs_src <- file.path(bin_path, pregs_name)

  if (!file.exists(pregs_src))
    stop("preGSf90 binary not found at: ", pregs_src,
         "\nInstall with install_genomic_programs().", call. = FALSE)

  # Copy binary to working directory to avoid DLL conflicts.
  # On Windows, the UGA genomic binaries are statically linked but will
  # pick up an incompatible libiomp5md.dll from the breedR bin directory
  # if run from there.
  pregs_bin <- file.path(dir, pregs_name)
  file.copy(pregs_src, pregs_bin, overwrite = TRUE)

  # Save and restore working directory
  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(pregs_bin)  # clean up copied binary
  })

  out <- system2(file.path(".", pregs_name), input = 'parameters',
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    # Save output for debugging before stopping
    writeLines(out, file.path(dir, "pregsf90.out"))
    stop("preGSf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 30), collapse = "\n"),
         call. = FALSE)
  }

  return(out)
}


#' Parse PREGSF90 QC output files
#'
#' Reads the quality control reports produced by PREGSF90 and returns
#' them as a structured R list.
#'
#' @param dir directory containing PREGSF90 output files.
#' @return A list with components: freq, excluded_snp, excluded_animals,
#'   conflicts, freq_raw, n_snp_total, n_snp_passed, n_snp_excluded,
#'   n_animals_excluded.
parse_pregsf90_qc <- function(dir) {

  result <- list()

  # Allele frequencies after QC
  freq_file <- file.path(dir, "freqdata.count.after.clean")
  if (file.exists(freq_file)) {
    freq <- utils::read.table(freq_file)
    names(freq) <- c("snp", "frequency", "exclusion_code")

    exclusion_labels <- c("Call Rate", "MAF", "Monomorphic",
                          "Excluded by request", "Mendelian error",
                          "HWE", "High correlation")
    freq$exclusion_reason <- ifelse(
      freq$exclusion_code == 0, "passed",
      exclusion_labels[freq$exclusion_code])

    result$freq <- freq
    result$excluded_snp <- freq[freq$exclusion_code > 0, ]
  }

  # Animals excluded by call rate
  callrate_file <- file.path(dir, "Gen_call_rate")
  if (file.exists(callrate_file) && file.info(callrate_file)$size > 0) {
    result$excluded_animals <- utils::read.table(callrate_file)
  }

  # Mendelian conflicts
  conflict_file <- file.path(dir, "Gen_conflicts")
  if (file.exists(conflict_file) && file.info(conflict_file)$size > 0) {
    result$conflicts <- readLines(conflict_file)
  }

  # Raw allele frequencies (before QC)
  raw_freq_file <- file.path(dir, "freqdata.count")
  if (file.exists(raw_freq_file)) {
    result$freq_raw <- utils::read.table(raw_freq_file)
    names(result$freq_raw) <- c("snp", "frequency")
  }

  # Summary statistics
  result$n_snp_total <- if (!is.null(result$freq)) nrow(result$freq) else NA
  result$n_snp_passed <- if (!is.null(result$freq))
    sum(result$freq$exclusion_code == 0) else NA
  result$n_snp_excluded <- if (!is.null(result$excluded_snp))
    nrow(result$excluded_snp) else 0L
  result$n_animals_excluded <- if (!is.null(result$excluded_animals))
    nrow(result$excluded_animals) else 0L

  return(result)
}


#' Format and write genotype data for BLUPF90
#'
#' Converts a genotype matrix from R into the fixed-width format required
#' by BLUPF90 programs (PREGSF90, BLUPF90+, RENUMF90). The format has two
#' fields per line: a right-justified animal ID and a concatenated genotype
#' string, with the genotype string starting at the same column on every line.
#'
#' Supports both integer genotypes (0/1/2 with 5 = missing) and fractional
#' genotypes from imputation (e.g., 0.50, 1.12). Fractional genotypes are
#' formatted with two decimal places and no separators (e.g., "0.501.120.25").
#'
#' @param geno numeric matrix. Rows = animals, columns = SNPs.
#'   For integer genotypes: values 0, 1, 2, or 5 (missing).
#'   For fractional genotypes: real values, typically from imputation.
#' @param ids character or integer vector. Animal IDs matching the pedigree.
#'   Must have the same length as \code{nrow(geno)}.
#' @param file character. Output file path.
#' @param fractional logical. If TRUE, write genotypes as real values with
#'   two decimal places (e.g., "1.000.502.00"). If FALSE (default), write
#'   as single digits (e.g., "10520").
#' @param missing_value integer. Code for missing genotypes (default 5).
#'   Only used when \code{fractional = FALSE}.
#'
#' @details
#' BLUPF90 requires that:
#' \enumerate{
#'   \item Fields are separated by at least one space
#'   \item The genotype string starts at the same column on every line
#'   \item For integer genotypes: single digits 0, 1, 2, 5 concatenated
#'   \item For fractional genotypes: values with exactly 2 decimal places,
#'     concatenated without separators
#' }
#'
#' @examples
#' \dontrun{
#' # Integer genotypes (SNP chip data)
#' geno <- matrix(c(0,1,2,1,0,2), nrow = 2)
#' write_snp_file(geno, ids = c("animal1", "animal2"), file = "snp.txt")
#' # Produces:
#' # animal1 012
#' # animal2 102
#'
#' # Fractional genotypes (imputed data)
#' geno_imp <- matrix(c(0.5, 1.12, 0.25, 1.5, 2.0, 0.0), nrow = 2)
#' write_snp_file(geno_imp, c("A1", "A2"), "snp_imp.txt", fractional = TRUE)
#' # Produces:
#' # A1 0.501.120.25
#' # A2 1.502.000.00
#' }
#' @export
write_snp_file <- function(geno, ids, file,
                            fractional = FALSE,
                            missing_value = 5L) {

  if (!is.matrix(geno)) geno <- as.matrix(geno)
  if (nrow(geno) != length(ids))
    stop("Number of rows in 'geno' (", nrow(geno),
         ") must match length of 'ids' (", length(ids), ").",
         call. = FALSE)

  ids_char <- as.character(ids)

  if (fractional) {
    # Fractional format: each value as X.XX (4 chars), no separators
    # e.g., 0.501.120.252.001.00
    geno_strings <- apply(geno, 1, function(row) {
      paste(sprintf("%.2f", row), collapse = "")
    })
  } else {
    # Integer format: single digit per SNP, concatenated
    # Replace NA with missing_value
    geno[is.na(geno)] <- missing_value
    geno_int <- matrix(as.integer(geno), nrow = nrow(geno))
    # Validate values
    valid <- c(0L, 1L, 2L, as.integer(missing_value))
    if (!all(geno_int %in% valid))
      warning("Genotype matrix contains values other than 0, 1, 2, ",
              missing_value, ". These will be written as-is.", call. = FALSE)
    geno_strings <- apply(geno_int, 1, paste, collapse = "")
  }

  # Right-justify IDs so genotypes start at a fixed column.
  # The ID field width is max(ID length) + 1 to ensure at least one
  # leading space even for the longest ID.
  id_width <- max(nchar(ids_char)) + 1L
  lines <- sprintf(paste0("%", id_width, "s %s"), ids_char, geno_strings)

  writeLines(lines, con = file)
  invisible(file)
}


#' Read a BLUPF90-format genotype file into R
#'
#' Reads a genotype file in BLUPF90 fixed-width format and returns the
#' IDs and genotype matrix as R objects.
#'
#' @param file character. Path to the genotype file.
#' @param fractional logical. If TRUE, parse genotypes as real values
#'   (4 chars each: X.XX). If FALSE (default), parse as single-digit integers.
#' @return A list with:
#'   \describe{
#'     \item{ids}{character vector of animal IDs}
#'     \item{geno}{numeric matrix of genotypes (rows = animals, cols = SNPs)}
#'   }
#' @export
read_snp_file <- function(file, fractional = FALSE) {
  raw <- readLines(file)
  raw <- raw[nchar(trimws(raw)) > 0]  # skip empty lines

  # Split each line into ID and genotype string
  # ID is everything before the genotype block (which starts at a fixed column)
  # Find the genotype start column from the first line
  first <- raw[1]
  # The genotype starts after the ID + space(s)
  # Detect: find last space before the genotype digits
  geno_start <- regexpr("[0-9]", sub("^\\s*\\S+\\s+", "", first))
  id_part <- trimws(sub("\\s+[0-9].*$", "", first))

  ids <- character(length(raw))
  geno_strings <- character(length(raw))

  for (i in seq_along(raw)) {
    parts <- strsplit(trimws(raw[i]), "\\s+", perl = TRUE)[[1]]
    ids[i] <- parts[1]
    geno_strings[i] <- paste(parts[-1], collapse = "")
  }

  if (fractional) {
    # Each value is 4 chars (X.XX)
    n_snp <- nchar(geno_strings[1]) / 4L
    geno <- matrix(NA_real_, nrow = length(ids), ncol = n_snp)
    for (i in seq_along(geno_strings)) {
      chars <- geno_strings[i]
      vals <- vapply(seq_len(n_snp), function(j) {
        as.numeric(substr(chars, (j - 1) * 4 + 1, j * 4))
      }, numeric(1))
      geno[i, ] <- vals
    }
  } else {
    # Each value is 1 char (single digit)
    n_snp <- nchar(geno_strings[1])
    geno <- matrix(NA_integer_, nrow = length(ids), ncol = n_snp)
    for (i in seq_along(geno_strings)) {
      geno[i, ] <- as.integer(strsplit(geno_strings[i], "")[[1]])
    }
  }

  return(list(ids = ids, geno = geno))
}


#' Write cross-reference ID file for PREGSF90
#'
#' Creates the XrefID file that maps renumbered pedigree IDs to original IDs.
#'
#' @param renumbered_ids integer vector. Renumbered IDs (from build_pedigree).
#' @param original_ids character or integer vector. Original animal IDs.
#' @param snp_file character. Path to the SNP genotype file. The XrefID file
#'   will be named \code{<snp_file>_XrefID}.
write_xref_file <- function(renumbered_ids, original_ids, snp_file) {
  xref <- data.frame(renumbered = renumbered_ids, original = original_ids)
  utils::write.table(xref, file = paste0(snp_file, "_XrefID"),
                      row.names = FALSE, col.names = FALSE, quote = FALSE)
}


### PostGSF90 — Genome-wide association study ###

#' Run PostGSF90 for GWAS analysis
#'
#' Extracts SNP effects and runs genome-wide association analysis on a
#' fitted genomic model. Must be called in the same R session as the
#' \code{\link{remlf90}} call that produced the model, since it reads
#' files from \code{tempdir()}.
#'
#' @param model a fitted \code{remlf90} object with a \code{genomic}
#'   component (i.e., fitted with \code{genomic = list(...)}).
#' @param windows_variance integer. Number of adjacent SNPs per window for
#'   computing variance explained. If NULL (default), not computed.
#' @param windows_variance_mbp numeric. Window size in megabases. If NULL
#'   (default), not computed.
#' @param windows_variance_type integer. 1 = moving windows (default),
#'   2 = exclusive (non-overlapping) windows.
#' @param which_weight integer or character. SNP weighting method for
#'   weighted ssGBLUP: 1 = y^2*2p(1-p), 2 = y^2, 4 or "nonlinearA" =
#'   VanRaden (2009). If NULL (default), not computed.
#' @param manhattan_plot logical. Generate Manhattan plot data (default TRUE).
#' @param snp_p_value logical. Compute p-values. Requires large memory
#'   (default FALSE).
#' @param postgs_trt_eff integer vector of length 2. Restrict analysis to
#'   specific trait and effect (c(trait, effect)). If NULL (default), all
#'   traits/effects.
#' @param extra_options character vector. Additional OPTION lines passed
#'   directly to PostGSF90.
#' @return A list with components:
#'   \describe{
#'     \item{snp_sol}{data.frame of SNP solutions, weights, and variance}
#'     \item{manhattan}{data.frame of Manhattan plot data (trait, effect,
#'       snp_effect, snp, chr, pos)}
#'     \item{pvalues}{data.frame of p-values (if snp_p_value = TRUE)}
#'     \item{windows}{data.frame of variance explained by windows}
#'     \item{snp_variance}{data.frame of per-SNP variance by chromosome}
#'     \item{output}{character vector of PostGSF90 stdout}
#'   }
#' @export
postgsf90 <- function(model,
                      windows_variance = NULL,
                      windows_variance_mbp = NULL,
                      windows_variance_type = 1L,
                      which_weight = NULL,
                      manhattan_plot = TRUE,
                      snp_p_value = FALSE,
                      postgs_trt_eff = NULL,
                      extra_options = NULL) {

  if (!inherits(model, 'remlf90'))
    stop("'model' must be a fitted remlf90 object.", call. = FALSE)
  if (is.null(model$genomic))
    stop("'model' must have been fitted with genomic = list(...).",
         call. = FALSE)

  # Files should be in tempdir() from the remlf90() call
  tmpdir <- tempdir()
  if (!file.exists(file.path(tmpdir, "solutions")))
    stop("Solutions file not found in tempdir(). ",
         "postgsf90() must be run in the same R session as remlf90().",
         call. = FALSE)

  # Build PostGSF90 OPTION lines
  postgs_opts <- build_postgsf90_options(
    windows_variance, windows_variance_mbp, windows_variance_type,
    which_weight, manhattan_plot, snp_p_value, postgs_trt_eff,
    extra_options
  )

  # Read the existing parameter file and rebuild with GWAS options.
  # Keep the base model specification, strip all OPTION lines,
  # then re-add only SNP_file, map_file, and GWAS options.
  par_lines <- readLines(file.path(tmpdir, "parameters"))
  base_lines <- par_lines[!grepl("^OPTION ", par_lines)]
  snp_opt <- par_lines[grepl("^OPTION SNP_file", par_lines)]
  map_opt <- par_lines[grepl("^OPTION map_file", par_lines)]
  new_par <- c(base_lines, snp_opt, map_opt,
               paste("OPTION", postgs_opts))
  writeLines(new_par, file.path(tmpdir, "parameters"))

  # Run PostGSF90
  bin_path <- breedR.getOption('breedR.bin')
  if (!check_genomic_programs(bin_path, quiet = TRUE))
    stop("Genomic program binaries (postGSf90) not installed. ",
         "See ?install_genomic_programs", call. = FALSE)

  postgs_out <- run_postgsf90(tmpdir, bin_path)
  writeLines(postgs_out, file.path(tmpdir, "postgsf90.out"))

  # Parse output files
  result <- parse_postgsf90(tmpdir)
  result$output <- postgs_out

  return(result)
}


## Build OPTION lines for PostGSF90
build_postgsf90_options <- function(windows_variance, windows_variance_mbp,
                                     windows_variance_type, which_weight,
                                     manhattan_plot, snp_p_value,
                                     postgs_trt_eff, extra_options) {
  opts <- character(0)

  if (!is.null(windows_variance))
    opts <- c(opts, paste("windows_variance", windows_variance))
  if (!is.null(windows_variance_mbp))
    opts <- c(opts, paste("windows_variance_mbp", windows_variance_mbp))
  if (!is.null(windows_variance) || !is.null(windows_variance_mbp))
    opts <- c(opts, paste("windows_variance_type", windows_variance_type))
  if (!is.null(which_weight))
    opts <- c(opts, paste("which_weight", which_weight))
  if (isTRUE(manhattan_plot))
    opts <- c(opts, "Manhattan_plot_R")
  if (isTRUE(snp_p_value))
    opts <- c(opts, "snp_p_value")
  if (!is.null(postgs_trt_eff)) {
    stopifnot(length(postgs_trt_eff) == 2)
    opts <- c(opts, paste("postgs_trt_eff",
                           postgs_trt_eff[1], postgs_trt_eff[2]))
  }
  if (!is.null(extra_options))
    opts <- c(opts, extra_options)

  return(opts)
}


#' Run PostGSF90 binary
#'
#' Executes the postGSf90 program in the specified directory.
#'
#' @param dir working directory containing parameter file, solutions,
#'   genotype file, and map file.
#' @param bin_path directory containing the postGSf90 binary.
#' @return Character vector of program stdout.
run_postgsf90 <- function(dir, bin_path) {

  postgs_name <- genomic_program_files(breedR.os.type())[2]
  postgs_src <- file.path(bin_path, postgs_name)

  if (!file.exists(postgs_src))
    stop("postGSf90 binary not found at: ", postgs_src,
         "\nInstall with install_genomic_programs().", call. = FALSE)

  # Copy binary to working directory (DLL isolation)
  postgs_bin <- file.path(dir, postgs_name)
  file.copy(postgs_src, postgs_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(postgs_bin)
  })

  out <- system2(file.path(".", postgs_name), input = 'parameters',
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "postgsf90.out"))
    stop("postGSf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 30), collapse = "\n"),
         call. = FALSE)
  }

  return(out)
}


#' Parse PostGSF90 output files
#'
#' Reads the GWAS results produced by PostGSF90 and returns them as
#' structured R data.frames.
#'
#' @param dir directory containing PostGSF90 output files.
#' @return A list with components: snp_sol, manhattan, pvalues, windows,
#'   window_segments, snp_variance.
parse_postgsf90 <- function(dir) {

  result <- list()

  # SNP solutions — always produced
  snp_sol_file <- file.path(dir, "snp_sol")
  if (file.exists(snp_sol_file) && file.info(snp_sol_file)$size > 0) {
    sol <- utils::read.table(snp_sol_file)
    base_names <- c("trait", "effect", "snp", "chr", "pos",
                    "solution", "weight")
    if (ncol(sol) >= 8) base_names <- c(base_names, "variance")
    if (ncol(sol) >= 9) base_names <- c(base_names, "var_solution")
    names(sol) <- base_names[seq_len(ncol(sol))]
    result$snp_sol <- sol
  }

  # Manhattan plot data
  chrsnp_file <- file.path(dir, "chrsnp")
  if (file.exists(chrsnp_file) && file.info(chrsnp_file)$size > 0) {
    mh <- utils::read.table(chrsnp_file)
    names(mh) <- c("trait", "effect", "snp_effect",
                    "snp", "chr", "pos")[seq_len(ncol(mh))]
    result$manhattan <- mh
  }

  # P-values
  pval_file <- file.path(dir, "chrsnp_pval")
  if (file.exists(pval_file) && file.info(pval_file)$size > 0) {
    pv <- utils::read.table(pval_file)
    names(pv) <- c("trait", "effect", "neg_log10_pval",
                    "snp", "chr", "pos")[seq_len(ncol(pv))]
    result$pvalues <- pv
  }

  # Variance by windows (non-overlapping, sums to 100%)
  winvar_file <- file.path(dir, "windows_variance")
  if (file.exists(winvar_file) && file.info(winvar_file)$size > 0) {
    wv <- utils::read.table(winvar_file, fill = TRUE)
    wv_names <- c("trait", "effect", "start_snp", "end_snp",
                   "window_size", "start_pos", "end_pos",
                   "window_id", "variance")
    names(wv) <- wv_names[seq_len(ncol(wv))]
    result$windows <- wv
  }

  # Window segment definitions
  winseg_file <- file.path(dir, "windows_segment")
  if (file.exists(winseg_file) && file.info(winseg_file)$size > 0) {
    ws <- utils::read.table(winseg_file, fill = TRUE)
    ws_names <- c("label", "window_size", "start_snp", "end_snp",
                   "window_id", "start_pos", "end_pos")
    names(ws) <- ws_names[seq_len(ncol(ws))]
    result$window_segments <- ws
  }

  # Per-SNP variance by chromosome (for plotting)
  chrsnpvar_file <- file.path(dir, "chrsnpvar")
  if (file.exists(chrsnpvar_file) && file.info(chrsnpvar_file)$size > 0) {
    cv <- utils::read.table(chrsnpvar_file)
    names(cv) <- c("trait", "effect", "variance",
                    "snp", "chr", "pos")[seq_len(ncol(cv))]
    result$snp_variance <- cv
  }

  # SNP predictions (allele freq + effects)
  snp_pred_file <- file.path(dir, "snp_pred")
  if (file.exists(snp_pred_file) && file.info(snp_pred_file)$size > 0) {
    result$snp_pred <- readLines(snp_pred_file)
  }

  return(result)
}


### PREDF90 — Predict DGV for new genotyped animals ###

#' Predict direct genomic values for new animals
#'
#' Uses SNP effects from \code{\link{postgsf90}} to predict direct genomic
#' values (DGV) for animals not in the original evaluation. This enables
#' genomic prediction for young/new animals based only on their genotypes.
#'
#' The prediction is: DGV = mu_hat + Z * a_hat, where a_hat are SNP effects
#' from PostGSF90 and Z is the centered genotype matrix.
#'
#' PREDF90 must be run in the same directory where PostGSF90 was run, since
#' it reads the \code{snp_pred} file produced by PostGSF90.
#'
#' @param snp_file character. Path to genotype file for animals to predict,
#'   in BLUPF90 format (see \code{\link{write_snp_file}}).
#' @param use_mu_hat logical. Add the base (mu_hat) to DGV so values are
#'   comparable to GEBV (default TRUE).
#' @param acc logical. Compute reliability of predictions (default FALSE).
#'   Requires \code{OPTION snp_p_value} in BLUPF90+ and
#'   \code{OPTION snp_var} in PostGSF90.
#' @param acc_type numeric. 1.0 for dairy cattle (reliability) or 0.5 for
#'   beef cattle (BIF accuracy). Default 1.0.
#' @param use_diagG_acc logical. Use inbreeding from G in the reliability
#'   denominator (default FALSE).
#' @param outfile character. Name of the output file (default
#'   "SNP_predictions").
#' @param dir character. Working directory containing PostGSF90 output files
#'   (snp_pred). Default \code{tempdir()}.
#' @return A data.frame with columns: id, call_rate, dgv, and optionally
#'   reliability.
#' @export
predf90 <- function(snp_file,
                    use_mu_hat = TRUE,
                    acc = FALSE,
                    acc_type = 1.0,
                    use_diagG_acc = FALSE,
                    outfile = "SNP_predictions",
                    dir = tempdir()) {

  if (!file.exists(snp_file))
    stop("SNP file not found: ", snp_file, call. = FALSE)

  # Check snp_pred file exists (from PostGSF90)
  if (!file.exists(file.path(dir, "snp_pred")))
    stop("snp_pred file not found in ", dir, ". ",
         "Run postgsf90() first.", call. = FALSE)

  bin_path <- breedR.getOption('breedR.bin')
  predf90_name <- "predf90"
  if (breedR.os.type() == 'windows') predf90_name <- "predf90.exe"

  predf90_src <- file.path(bin_path, predf90_name)
  if (!file.exists(predf90_src))
    stop("predf90 binary not found. Run install_genomic_programs().",
         call. = FALSE)

  # Copy SNP file to working directory if not already there
  snp_basename <- basename(snp_file)
  snp_dest <- file.path(dir, snp_basename)
  if (normalizePath(snp_file, mustWork = FALSE) !=
      normalizePath(snp_dest, mustWork = FALSE))
    file.copy(snp_file, snp_dest, overwrite = TRUE)

  # Copy binary to working directory (DLL isolation)
  predf90_bin <- file.path(dir, predf90_name)
  file.copy(predf90_src, predf90_bin, overwrite = TRUE)

  # Build command-line arguments
  args <- c("--snpfile", snp_basename)
  if (isTRUE(use_mu_hat)) args <- c(args, "--use_mu_hat")
  if (isTRUE(acc)) {
    args <- c(args, "--acc")
    args <- c(args, "--acc_type", as.character(acc_type))
  }
  if (isTRUE(use_diagG_acc)) args <- c(args, "--use_diagG_acc")
  if (!is.null(outfile)) args <- c(args, "--outfile", outfile)

  # Execute
  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(predf90_bin)
  })

  out <- system2(file.path(".", predf90_name), args = args,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "predf90.out"))
    stop("predf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(out, file.path(dir, "predf90.out"))

  # Parse output
  out_path <- file.path(dir, outfile)
  if (!file.exists(out_path))
    stop("predf90 did not produce output file: ", outfile, call. = FALSE)

  predictions <- utils::read.table(out_path, header = FALSE)
  col_names <- c("id", "call_rate", "dgv")
  if (ncol(predictions) >= 4) col_names <- c(col_names, "reliability")
  names(predictions) <- col_names[seq_len(ncol(predictions))]

  return(predictions)
}
