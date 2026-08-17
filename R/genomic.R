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

  # Save the plain genomic G-inverse so a later postgsf90() GWAS run can read
  # it (postGSf90 needs the G-inverse via readGInverse, not the ssGBLUP
  # correction matrix GimA22i). Opt-in, since the file is large.
  if (is.null(genomic$save_ginverse)) genomic$save_ginverse <- FALSE
  if (!is.logical(genomic$save_ginverse) || length(genomic$save_ginverse) != 1)
    stop("'genomic$save_ginverse' must be a single logical.", call. = FALSE)

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
  # Save the plain G-inverse for a downstream postgsf90() GWAS run.
  if (isTRUE(genomic$save_ginverse)) opts <- c(opts, "saveGInverse")

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
#' @param pedigree the renumbered pedigree (from the fitted genetic effect),
#'   used to derive the SNP cross-reference file. May be NULL.
#' @return List with updated \code{pf90} and \code{pregs_out}.
run_pregsf90_pipeline <- function(genomic, genomic_opts, pf90,
                                   tmpdir, breedR.bin, pedigree = NULL) {
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

  # preGSf90 needs a cross-reference file mapping each SNP-file (original)
  # animal ID to the renumbered pedigree code that breedR writes in the
  # pedigree/data files. Use a user-supplied XrefID if present; otherwise
  # derive it from the pedigree's renumbering.
  xref_file <- paste0(genomic$snp_file, "_XrefID")
  if (file.exists(xref_file)) {
    copy_if_needed(xref_file)
  } else if (!is.null(pedigree)) {
    write_xref_from_pedigree(genomic$snp_file, pedigree, tmpdir)
  } else {
    stop("Cannot build the genomic XrefID: no pedigree available and no ",
         "'", basename(xref_file), "' provided.", call. = FALSE)
  }

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

  # Allele frequencies after QC. preGSf90 writes two columns (snp index,
  # frequency); some versions/options add a third exclusion-code column. Handle
  # both rather than assume three (which errors when only two are present).
  freq_file <- file.path(dir, "freqdata.count.after.clean")
  if (file.exists(freq_file)) {
    freq <- utils::read.table(freq_file)
    if (ncol(freq) >= 3) {
      names(freq)[1:3] <- c("snp", "frequency", "exclusion_code")

      # Map exclusion codes to labels by name (code 0 = passed). Indexing a
      # plain vector with code 0 silently drops that element and misaligns
      # every subsequent label, so key the lookup by the code itself.
      exclusion_labels <- c("0" = "passed",
                            "1" = "Call Rate", "2" = "MAF", "3" = "Monomorphic",
                            "4" = "Excluded by request", "5" = "Mendelian error",
                            "6" = "HWE", "7" = "High correlation")
      freq$exclusion_reason <- unname(
        exclusion_labels[as.character(freq$exclusion_code)])
      result$excluded_snp <- freq[freq$exclusion_code > 0, ]
    } else {
      names(freq)[1:2] <- c("snp", "frequency")
    }
    result$freq <- freq
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

  # Summary statistics. freqdata.count.after.clean lists the SNPs that passed
  # QC; freqdata.count is the raw set. With the 2-column format there is no
  # per-SNP exclusion code, so passed = rows in the after-clean file and the
  # totals come from the raw file.
  n_raw  <- if (!is.null(result$freq_raw)) nrow(result$freq_raw) else NULL
  n_pass <- if (!is.null(result$freq)) {
    if (!is.null(result$freq$exclusion_code))
      sum(result$freq$exclusion_code == 0)
    else nrow(result$freq)
  } else NULL

  result$n_snp_total <- if (!is.null(n_raw)) n_raw
    else if (!is.null(n_pass)) n_pass else NA
  result$n_snp_passed <- if (!is.null(n_pass)) n_pass else NA
  result$n_snp_excluded <- if (!is.null(result$excluded_snp))
      nrow(result$excluded_snp)
    else if (!is.null(n_raw) && !is.null(n_pass)) n_raw - n_pass
    else 0L
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
    # Fractional format: each value as X.XX (exactly 4 chars), no separators
    # e.g., 0.501.120.252.001.00
    # The reader parses this as fixed 4-char fields, so anything that does not
    # format to 4 characters (NA -> "NA", negatives, values >= 10) would corrupt
    # the layout silently. Reject those inputs rather than write an unparseable
    # file. (A dedicated missing-value convention is pending binary verification.)
    if (anyNA(geno))
      stop("Fractional SNP format does not support missing (NA) genotypes.",
           call. = FALSE)
    geno_strings <- apply(geno, 1, function(row) {
      tok <- sprintf("%.2f", row)
      bad <- which(nchar(tok) != 4L)
      if (length(bad))
        stop("Fractional genotype values must format to exactly 4 characters ",
             "(0.00-9.99); got '", tok[bad[1]], "'.", call. = FALSE)
      paste(tok, collapse = "")
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

  ids <- character(length(raw))
  geno_strings <- character(length(raw))

  for (i in seq_along(raw)) {
    parts <- strsplit(trimws(raw[i]), "\\s+", perl = TRUE)[[1]]
    ids[i] <- parts[1]
    geno_strings[i] <- paste(parts[-1], collapse = "")
  }

  # All genotype rows must share the same length, otherwise the fixed-width
  # parsing below would silently return a mis-shaped matrix.
  widths <- nchar(geno_strings)
  if (length(unique(widths)) != 1L)
    stop("Genotype rows have inconsistent lengths (", min(widths), "-",
         max(widths), " chars); cannot parse SNP file.", call. = FALSE)

  if (fractional) {
    # Each value is 4 chars (X.XX)
    if (widths[1] %% 4L != 0L)
      stop("Fractional SNP line length (", widths[1],
           ") is not a multiple of 4.", call. = FALSE)
    n_snp <- widths[1] / 4L
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


## Derive the preGSf90 XrefID from breedR's pedigree renumbering and write it
## into the working directory next to the copied SNP file. Each SNP-file row's
## (original) animal id is mapped to its renumbered code (its position in the
## pedigree labels). Errors if a genotyped animal is absent from the pedigree.
write_xref_from_pedigree <- function(snp_file, pedigree, tmpdir) {
  raw <- readLines(snp_file)
  raw <- raw[nchar(trimws(raw)) > 0]
  snp_ids <- vapply(strsplit(trimws(raw), "[[:space:]]+"), `[`, "", 1L)

  labels <- as.character(pedigree@label)
  renum <- match(snp_ids, labels)
  if (anyNA(renum))
    stop("Genotyped animals not found in the pedigree: ",
         paste(utils::head(snp_ids[is.na(renum)], 10L), collapse = ", "),
         if (sum(is.na(renum)) > 10L) ", ..." else "", call. = FALSE)

  utils::write.table(
    data.frame(renum, snp_ids),
    file = file.path(tmpdir, paste0(basename(snp_file), "_XrefID")),
    row.names = FALSE, col.names = FALSE, quote = FALSE)
}


### PostGSF90 — Genome-wide association study ###

#' Run PostGSF90 for GWAS analysis
#'
#' Back-solves the genomic breeding values into SNP effects and runs a
#' genome-wide association analysis on a fitted single-step model. Must be
#' called in the same R session as the \code{\link{remlf90}} call that produced
#' the model, since it reads that fit's working directory
#' (\code{model$reml$dir}), which lives under \code{tempdir()}.
#'
#' The model must have been fitted with \code{save_ginverse = TRUE} in its
#' \code{genomic} list; postGSf90 reads the saved genomic G-inverse. For
#' chromosome/position information in the Manhattan output, supply a
#' \code{map_file} at fit time.
#'
#' @param model a fitted \code{remlf90} object fitted with
#'   \code{genomic = list(..., save_ginverse = TRUE)}.
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
#'     \item{dir}{the working directory the SNP effects were written to. Pass
#'       it to \code{\link{predf90}} as its \code{dir}.}
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

  # postGSf90 back-solves SNP effects from the plain genomic G-inverse, which is
  # only written when the fit was run with save_ginverse = TRUE. Without it the
  # working directory has no G-inverse to read.
  if (!isTRUE(model$genomic$save_ginverse))
    stop("This model was not fitted for GWAS. Refit with ",
         "genomic = list(..., save_ginverse = TRUE) before calling postgsf90().",
         call. = FALSE)

  ## Read the working directory of *this* model rather than the session's.
  ## Objects rebuilt by breedR.qget() predate the recording, hence the
  ## fallback; a fit from this session always carries its own path.
  tmpdir <- if (!is.null(model$reml$dir)) model$reml$dir else tempdir()
  if (!file.exists(file.path(tmpdir, "solutions")))
    stop("Solutions file not found in ", tmpdir, ". ",
         "postgsf90() must be run in the same R session as remlf90().",
         call. = FALSE)

  # Build PostGSF90 OPTION lines
  postgs_opts <- build_postgsf90_options(
    windows_variance, windows_variance_mbp, windows_variance_type,
    which_weight, manhattan_plot, snp_p_value, postgs_trt_eff,
    extra_options
  )

  # Rebuild the parameter file for GWAS. postGSf90 reads the plain G-inverse
  # (readGInverse) that the fit saved via saveGInverse; it does NOT use the
  # blupf90+ ssGBLUP option readGimA22i (that is the G-inverse minus A22-inverse
  # correction, a different matrix). Keep SNP_file, translate the map to
  # chrinfo, drop the REML-solver options and readGimA22i, then add the G-inverse
  # read option and the GWAS options.
  par_lines <- readLines(file.path(tmpdir, "parameters"))
  drop_opts <- paste0("^OPTION (sol se|method|EM-REML|se_covar_function|",
                      "maxrounds|conv_crit|use_yams|readGimA22i)\\b")
  kept_lines <- par_lines[!grepl(drop_opts, par_lines)]
  # map_file -> chrinfo (postGSf90 uses chrinfo for chromosome/position info)
  kept_lines <- sub("^OPTION map_file ", "OPTION chrinfo ", kept_lines)
  new_par <- c(kept_lines, "OPTION readGInverse", paste("OPTION", postgs_opts))
  ## Write it beside the fit's parameter file rather than over it. That file is
  ## part of what the model object now advertises through res$reml$dir, and
  ## rewriting it in place both destroyed it and made a second postgsf90() on
  ## the same model wrong -- drop_opts above does not match readGInverse or the
  ## GWAS options, so they accumulated on every call.
  postgs_par <- "parameters_postgs"
  writeLines(new_par, file.path(tmpdir, postgs_par))

  # postGSf90 reads external inbreeding coefficients (renf90.inb) when a
  # pedigree is present. breedR's ssGBLUP fit does not track inbreeding (A22 is
  # built with F = 0 via readGimA22i), so write a matching zero-inbreeding file
  # for the pedigree animals; otherwise postGSf90 aborts ("inbreeding file not
  # found").
  ped_file <- file.path(tmpdir, "pedigree_genetic")
  if (file.exists(ped_file)) {
    n_anim <- length(readLines(ped_file))
    utils::write.table(
      data.frame(id = seq_len(n_anim), f = 0),
      file.path(tmpdir, "renf90.inb"),
      row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  # Run PostGSF90
  bin_path <- breedR.getOption('breedR.bin')
  if (!check_genomic_programs(bin_path, quiet = TRUE))
    stop("Genomic program binaries (postGSf90) not installed. ",
         "See ?install_genomic_programs", call. = FALSE)

  postgs_out <- run_postgsf90(tmpdir, bin_path, par_name = postgs_par)
  writeLines(postgs_out, file.path(tmpdir, "postgsf90.out"))

  # postGSf90 can exit 0 (or write a stub snp_sol) while its log reports an
  # error, so verify the run actually succeeded rather than returning a
  # garbled/empty result.
  # "ERROR:" matches the Fortran fatal-error convention without false-matching
  # benign output like "number of errors = 0" or "standard error".
  snp_sol_file <- file.path(tmpdir, "snp_sol")
  if (any(grepl("ERROR:", postgs_out)) ||
      !file.exists(snp_sol_file) ||
      file.info(snp_sol_file)$size == 0)
    stop("postGSf90 failed to produce valid SNP solutions. Check the log:\n",
         file.path(tmpdir, "postgsf90.out"), call. = FALSE)

  # Parse output files
  result <- parse_postgsf90(tmpdir)
  result$output <- postgs_out
  ## Where the SNP effects landed, so predf90() can be pointed at them without
  ## guessing at the session's tempdir().
  result$dir <- tmpdir

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
#' @param par_name name of the parameter file to feed the program, relative to
#'   \code{dir}. Defaults to the GWAS parameter file postgsf90() writes, which
#'   is kept separate from the fit's own \code{parameters}.
#' @return Character vector of program stdout.
run_postgsf90 <- function(dir, bin_path, par_name = 'parameters_postgs') {

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

  out <- system2(file.path(".", postgs_name), input = par_name,
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

  # SNP solutions — always produced. postGSf90 writes a header row
  # (trait effect snp chr pos snp_effect weight variance_explained var_a_hat).
  snp_sol_file <- file.path(dir, "snp_sol")
  if (file.exists(snp_sol_file) && file.info(snp_sol_file)$size > 0) {
    sol <- utils::read.table(snp_sol_file, header = TRUE)
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
#' @param pedfile character or NULL. Pedigree file enabling reliability of
#'   prediction from the pedigree (RPG). When NULL (default), predf90 is run
#'   with \code{--no_rpg} and only genomic DGV are produced; predf90 aborts
#'   without writing predictions if neither a pedigree nor \code{--no_rpg} is
#'   supplied.
#' @param outfile character. Name of the output file (default
#'   "SNP_predictions").
#' @param dir character. Working directory containing PostGSF90 output files
#'   (snp_pred). Pass the \code{dir} element of the \code{\link{postgsf90}}
#'   result: each fit now works in its own subdirectory of \code{tempdir()},
#'   so the default is only right when no model was fitted in this session.
#' @return A data.frame with columns: id, call_rate, dgv, and optionally
#'   reliability.
#' @export
predf90 <- function(snp_file,
                    use_mu_hat = TRUE,
                    acc = FALSE,
                    acc_type = 1.0,
                    use_diagG_acc = FALSE,
                    pedfile = NULL,
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

  # predf90 needs either a pedigree file (for reliability from the pedigree,
  # RPG) or --no_rpg; without one it computes DGV but exits without writing the
  # output file. This wrapper predicts genomic DGV, so default to --no_rpg.
  if (!is.null(pedfile)) {
    if (!file.exists(pedfile))
      stop("Pedigree file not found: ", pedfile, call. = FALSE)
    ped_dest <- file.path(dir, basename(pedfile))
    if (normalizePath(pedfile, mustWork = FALSE) !=
        normalizePath(ped_dest, mustWork = FALSE))
      file.copy(pedfile, ped_dest, overwrite = TRUE)
    args <- c(args, "--pedfile", basename(pedfile))
  } else {
    args <- c(args, "--no_rpg")
  }

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

  if (any(grepl("ERROR:", out)))
    stop("predf90 reported an error despite a clean exit. Check the log:\n",
         file.path(dir, "predf90.out"), call. = FALSE)

  # Parse output
  out_path <- file.path(dir, outfile)
  if (!file.exists(out_path))
    stop("predf90 did not produce output file '", outfile,
         "' despite a clean exit. Last lines of its output:\n",
         paste(utils::tail(out, 8), collapse = "\n"),
         "\nFull log: ", file.path(dir, "predf90.out"), call. = FALSE)

  predictions <- utils::read.table(out_path, header = FALSE)
  col_names <- c("id", "call_rate", "dgv")
  if (ncol(predictions) >= 4) col_names <- c(col_names, "reliability")
  names(predictions) <- col_names[seq_len(ncol(predictions))]

  return(predictions)
}
