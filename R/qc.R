### QCF90 — Genotype and pedigree quality control ###


#' Quality control for genotype and pedigree data
#'
#' Runs QCF90 on genotype and optionally pedigree data to perform
#' quality checks including call rate filtering, MAF filtering,
#' monomorphic marker removal, Hardy-Weinberg equilibrium testing,
#' Mendelian conflict detection, and duplicate sample identification.
#'
#' QCF90 works with non-renumbered (raw alphanumeric) files, making it
#' suitable as a first-pass diagnostic before \code{\link{renumf90}}.
#'
#' @param snp_file character. Path to genotype file in BLUPF90 format.
#' @param ped_file character or NULL. Path to pedigree file
#'   (animal, sire, dam — space-separated).
#' @param map_file character or NULL. Path to SNP map file
#'   (SNP_ID, CHR, POS header).
#' @param call_rate_markers numeric. Minimum call rate for markers
#'   (default 0.90). Must be between 0 and 1.
#' @param call_rate_animals numeric. Minimum call rate for animals
#'   (default 0.90). Must be between 0 and 1.
#' @param maf numeric. Minimum allele frequency threshold (default 0.05).
#'   Must be between 0 and 0.5.
#' @param hwe numeric or NULL. HWE deviation threshold. NULL skips the
#'   check (default NULL).
#' @param check_parentage logical. Check Mendelian inconsistencies.
#'   Defaults to TRUE when \code{ped_file} is provided, FALSE otherwise.
#' @param trio logical. Use sire-dam-progeny trio for Mendelian checks
#'   (default FALSE).
#' @param check_identity logical or character. Check for duplicate samples.
#'   If character, path to a file of IDs to compare against (default FALSE).
#' @param identity_threshold numeric. Threshold for identity detection
#'   (default 0.99).
#' @param sex_chr integer vector or NULL. Chromosome IDs for sex chromosomes.
#' @param skip_chr integer vector or NULL. Chromosomes to skip in QC.
#' @param exclude_chr integer vector or NULL. Chromosomes to remove entirely.
#' @param save_clean logical. Save cleaned genotype files with unqualified
#'   markers/animals removed (default TRUE).
#' @param save_log logical. Save QC log file (default TRUE).
#' @param outcallrate logical. Output detailed call rate per marker and
#'   animal (default FALSE).
#' @param remove_markers logical. Remove unqualified markers progressively
#'   during QC steps (default TRUE).
#' @param remove_animals logical. Remove unqualified animals progressively
#'   (default TRUE).
#' @param extra_args character vector or NULL. Additional command-line
#'   arguments passed directly to QCF90. Use for options not covered by
#'   other parameters (e.g., \code{c("--no-check-format", "--fastread")}).
#' @param dir character. Working directory for output files.
#' @return A list with:
#'   \describe{
#'     \item{log}{character vector of QC log content}
#'     \item{status}{character vector of QC status file content}
#'     \item{clean_snp}{path to cleaned SNP file (if save_clean)}
#'     \item{clean_xref}{path to cleaned XrefID file (if present)}
#'     \item{removed_snp}{character vector of removed SNP info}
#'     \item{removed_animals}{character vector of removed animal info}
#'     \item{output}{character vector of program stdout}
#'   }
#' @export
qcf90 <- function(snp_file,
                   ped_file = NULL,
                   map_file = NULL,
                   call_rate_markers = 0.90,
                   call_rate_animals = 0.90,
                   maf = 0.05,
                   hwe = NULL,
                   check_parentage = NULL,
                   trio = FALSE,
                   check_identity = FALSE,
                   identity_threshold = 0.99,
                   sex_chr = NULL,
                   skip_chr = NULL,
                   exclude_chr = NULL,
                   save_clean = TRUE,
                   save_log = TRUE,
                   outcallrate = FALSE,
                   remove_markers = TRUE,
                   remove_animals = TRUE,
                   extra_args = NULL,
                   dir = tempdir()) {

  # Validate inputs
  if (!file.exists(snp_file))
    stop("SNP file not found: ", snp_file, call. = FALSE)
  if (!is.null(ped_file) && !file.exists(ped_file))
    stop("Pedigree file not found: ", ped_file, call. = FALSE)
  if (!is.null(map_file) && !file.exists(map_file))
    stop("Map file not found: ", map_file, call. = FALSE)

  # Validate numeric thresholds
  if (!is.numeric(call_rate_markers) || call_rate_markers < 0 ||
      call_rate_markers > 1)
    stop("'call_rate_markers' must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(call_rate_animals) || call_rate_animals < 0 ||
      call_rate_animals > 1)
    stop("'call_rate_animals' must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(maf) || maf < 0 || maf > 0.5)
    stop("'maf' must be between 0 and 0.5.", call. = FALSE)
  if (!is.null(hwe) && (!is.numeric(hwe) || hwe < 0 || hwe > 1))
    stop("'hwe' must be between 0 and 1, or NULL.", call. = FALSE)

  # Default check_parentage to TRUE when pedigree is provided
  if (is.null(check_parentage))
    check_parentage <- !is.null(ped_file)

  bin_path <- breedR.getOption('breedR.bin')
  qc_name <- "qcf90"
  if (breedR.os.type() == 'windows') qc_name <- paste0(qc_name, ".exe")
  qc_src <- file.path(bin_path, qc_name)

  if (!file.exists(qc_src))
    stop("qcf90 binary not found. Run install_genomic_programs().",
         call. = FALSE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Copy input files to working directory
  copy_if_needed <- function(src) {
    dest <- file.path(dir, basename(src))
    if (normalizePath(src, mustWork = FALSE) !=
        normalizePath(dest, mustWork = FALSE))
      file.copy(src, dest, overwrite = TRUE)
  }
  copy_if_needed(snp_file)
  if (!is.null(ped_file))
    copy_if_needed(ped_file)
  if (!is.null(map_file))
    copy_if_needed(map_file)

  # Copy identity comparison file if it's a path
  if (is.character(check_identity) && file.exists(check_identity))
    copy_if_needed(check_identity)

  # Build command-line arguments
  args <- build_qcf90_args(
    snp_file = snp_file, ped_file = ped_file, map_file = map_file,
    call_rate_markers = call_rate_markers,
    call_rate_animals = call_rate_animals,
    maf = maf, hwe = hwe, check_parentage = check_parentage,
    trio = trio, check_identity = check_identity,
    identity_threshold = identity_threshold,
    sex_chr = sex_chr, skip_chr = skip_chr, exclude_chr = exclude_chr,
    save_clean = save_clean, save_log = save_log,
    outcallrate = outcallrate, remove_markers = remove_markers,
    remove_animals = remove_animals, extra_args = extra_args
  )

  # Copy binary (DLL isolation)
  qc_bin <- file.path(dir, qc_name)
  file.copy(qc_src, qc_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(qc_bin)
  })

  out <- system2(file.path(".", qc_name), args = args,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "qcf90_stdout.txt"))
    stop("qcf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(out, file.path(dir, "qcf90_stdout.txt"))

  # Parse output files
  result <- list(output = out)

  log_file <- file.path(dir, "qcf90.log")
  if (file.exists(log_file))
    result$log <- readLines(log_file)

  status_file <- file.path(dir, "qcf90.status")
  if (file.exists(status_file))
    result$status <- readLines(status_file)

  # Clean files — search for any file matching the _clean pattern
  snp_basename <- basename(snp_file)
  clean_files <- list.files(dir, pattern = "_clean$", full.names = TRUE)
  if (length(clean_files) > 0)
    result$clean_snp <- clean_files[1]

  xref_clean_files <- list.files(dir, pattern = "_clean_XrefID$",
                                  full.names = TRUE)
  if (length(xref_clean_files) > 0)
    result$clean_xref <- xref_clean_files[1]

  # Removed files
  removed_snp_files <- list.files(dir, pattern = "_SNPs_removed$",
                                   full.names = TRUE)
  if (length(removed_snp_files) > 0)
    result$removed_snp <- readLines(removed_snp_files[1])

  removed_anim_files <- list.files(dir, pattern = "_Animals_removed$",
                                    full.names = TRUE)
  if (length(removed_anim_files) > 0)
    result$removed_animals <- readLines(removed_anim_files[1])

  return(result)
}


## Build QCF90 command-line arguments from R parameters.
## Each chromosome number is passed as a separate argument element
## so system2() tokenizes them correctly.
build_qcf90_args <- function(snp_file, ped_file, map_file,
                              call_rate_markers, call_rate_animals,
                              maf, hwe, check_parentage, trio,
                              check_identity, identity_threshold,
                              sex_chr, skip_chr, exclude_chr,
                              save_clean, save_log, outcallrate,
                              remove_markers, remove_animals,
                              extra_args) {

  args <- c("--snpfile", basename(snp_file))

  if (!is.null(ped_file))
    args <- c(args, "--pedfile", basename(ped_file))
  if (!is.null(map_file))
    args <- c(args, "--mapfile", basename(map_file))

  # QC thresholds
  args <- c(args, "--crm", as.character(call_rate_markers))
  args <- c(args, "--cra", as.character(call_rate_animals))
  args <- c(args, "--maf", as.character(maf))

  if (!is.null(hwe))
    args <- c(args, "--hwe", as.character(hwe))

  # Parentage
  if (isTRUE(check_parentage))
    args <- c(args, "--check-parentage")
  if (isTRUE(trio))
    args <- c(args, "--trio")

  # Identity — copy basename if it's a file path
  if (!isFALSE(check_identity)) {
    if (is.character(check_identity))
      args <- c(args, "--check-identity", basename(check_identity))
    else
      args <- c(args, "--check-identity")
    args <- c(args, "--threshold-identity", as.character(identity_threshold))
  }

  # Chromosome filtering — each number as a separate arg element
  if (!is.null(sex_chr))
    args <- c(args, "--sex-chr", as.character(sex_chr))
  if (!is.null(skip_chr))
    args <- c(args, "--skip-chr", as.character(skip_chr))
  if (!is.null(exclude_chr))
    args <- c(args, "--exclude-chr", as.character(exclude_chr))

  # Output options
  if (isTRUE(save_clean))
    args <- c(args, "--save-clean")
  if (isTRUE(save_log))
    args <- c(args, "--save-log")
  if (isTRUE(outcallrate))
    args <- c(args, "--outcallrate")
  if (isTRUE(remove_markers))
    args <- c(args, "--remove-markers")
  if (isTRUE(remove_animals))
    args <- c(args, "--remove-animals")

  if (!is.null(extra_args))
    args <- c(args, extra_args)

  return(args)
}
