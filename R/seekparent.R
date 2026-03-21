### SEEKPARENTF90 — Parentage verification and discovery ###


#' Verify and discover parentage using SNP markers
#'
#' Uses genotype data to check parent-offspring relationships for Mendelian
#' conflicts and optionally searches for the correct parent among a set of
#' candidates. Based on Hayes (2011) and Wiggans et al. (2010).
#'
#' Supports alphanumeric IDs, multiple SNP chips, and year-of-birth
#' filtering for candidate parents.
#'
#' @param snp_file character. Path to genotype file in BLUPF90 format
#'   (see \code{\link{write_snp_file}}).
#' @param ped_file character. Path to pedigree file (animal, sire, dam,
#'   space-separated). Use 0 for unknown parents.
#' @param yob logical. If TRUE, read year of birth from the 4th column of
#'   the pedigree file (default FALSE).
#' @param seek_sire logical or character. If TRUE, search for correct sire
#'   among genotyped sires in pedigree. If a file path, use that file as
#'   the candidate sire list (default FALSE).
#' @param seek_dam logical or character. If TRUE, search for correct dam
#'   among genotyped dams in pedigree. If a file path, use that file as
#'   the candidate dam list (default FALSE).
#' @param seek_type integer. Which animals to search: 1 = only non-matching
#'   parents (default), 2 = all genotyped individuals.
#' @param only_in_list character or NULL. Path to a file listing animal IDs
#'   to restrict the analysis to.
#' @param excl_thr_prob numeric. Exclusion threshold as percentage of SNPs
#'   (default 1.0 for chips > 130 SNPs).
#' @param assign_thr_prob numeric. Assignment threshold as percentage
#'   (default 0.5).
#' @param call_rate numeric. Minimum call rate to include a sample
#'   (default 0.90).
#' @param trio logical. Use both sire and dam info for Mendelian conflict
#'   detection, allowing non-homozygous markers (default FALSE).
#' @param chips character or NULL. Path to chip map file for multi-chip
#'   analyses.
#' @param check_duplicates logical or numeric. Check for duplicate samples.
#'   If numeric, use as the threshold (default FALSE, threshold 0.9).
#' @param extra_args character vector. Additional command-line arguments.
#' @param dir character. Working directory for output files.
#' @return A list with:
#'   \describe{
#'     \item{check}{data.frame of parent-offspring check results}
#'     \item{assigned}{data.frame of corrected pedigree (if seek was used)}
#'     \item{seek_sire}{data.frame of sire search results (if requested)}
#'     \item{seek_dam}{data.frame of dam search results (if requested)}
#'     \item{output}{character vector of program stdout}
#'   }
#' @export
seekparentf90 <- function(snp_file,
                           ped_file,
                           yob = FALSE,
                           seek_sire = FALSE,
                           seek_dam = FALSE,
                           seek_type = 1L,
                           only_in_list = NULL,
                           excl_thr_prob = NULL,
                           assign_thr_prob = NULL,
                           call_rate = NULL,
                           trio = FALSE,
                           chips = NULL,
                           check_duplicates = FALSE,
                           extra_args = NULL,
                           dir = tempdir()) {

  if (!file.exists(snp_file))
    stop("SNP file not found: ", snp_file, call. = FALSE)
  if (!file.exists(ped_file))
    stop("Pedigree file not found: ", ped_file, call. = FALSE)

  bin_path <- breedR.getOption('breedR.bin')
  seek_name <- "seekparentf90"
  if (breedR.os.type() == 'windows') seek_name <- paste0(seek_name, ".exe")
  seek_src <- file.path(bin_path, seek_name)

  if (!file.exists(seek_src))
    stop("seekparentf90 binary not found. Run install_genomic_programs().",
         call. = FALSE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Copy input files to working directory
  snp_dest <- file.path(dir, basename(snp_file))
  ped_dest <- file.path(dir, basename(ped_file))
  if (normalizePath(snp_file, mustWork = FALSE) !=
      normalizePath(snp_dest, mustWork = FALSE))
    file.copy(snp_file, snp_dest, overwrite = TRUE)
  if (normalizePath(ped_file, mustWork = FALSE) !=
      normalizePath(ped_dest, mustWork = FALSE))
    file.copy(ped_file, ped_dest, overwrite = TRUE)

  # Build command-line arguments
  args <- c("--pedfile", basename(ped_file),
            "--snpfile", basename(snp_file))

  if (isTRUE(yob)) args <- c(args, "--yob")

  # Sire seeking
  if (isTRUE(seek_sire)) {
    args <- c(args, "--seeksire_in_ped")
  } else if (is.character(seek_sire) && file.exists(seek_sire)) {
    file.copy(seek_sire, file.path(dir, basename(seek_sire)),
              overwrite = TRUE)
    args <- c(args, "--seeksire", basename(seek_sire))
  }

  # Dam seeking
  if (isTRUE(seek_dam)) {
    args <- c(args, "--seekdam_in_ped")
  } else if (is.character(seek_dam) && file.exists(seek_dam)) {
    file.copy(seek_dam, file.path(dir, basename(seek_dam)),
              overwrite = TRUE)
    args <- c(args, "--seekdam", basename(seek_dam))
  }

  if (seek_type != 1L) args <- c(args, "--seektype", as.character(seek_type))

  if (!is.null(only_in_list)) {
    if (file.exists(only_in_list)) {
      file.copy(only_in_list, file.path(dir, basename(only_in_list)),
                overwrite = TRUE)
      args <- c(args, "--only_in_list", basename(only_in_list))
    }
  }

  if (!is.null(excl_thr_prob))
    args <- c(args, "--excl_thr_prob", as.character(excl_thr_prob))
  if (!is.null(assign_thr_prob))
    args <- c(args, "--assign_thr_prob", as.character(assign_thr_prob))
  if (!is.null(call_rate))
    args <- c(args, "--thr_call_rate", as.character(call_rate))
  if (isTRUE(trio)) args <- c(args, "--trio")

  if (!is.null(chips)) {
    file.copy(chips, file.path(dir, basename(chips)), overwrite = TRUE)
    args <- c(args, "--chips", basename(chips))
  }

  if (!isFALSE(check_duplicates)) {
    if (is.numeric(check_duplicates))
      args <- c(args, "--duplicate", as.character(check_duplicates))
    else
      args <- c(args, "--duplicate")
  }

  if (!is.null(extra_args)) args <- c(args, extra_args)

  # Copy binary (DLL isolation)
  seek_bin <- file.path(dir, seek_name)
  file.copy(seek_src, seek_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(seek_bin)
  })

  out <- system2(file.path(".", seek_name), args = args,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "seekparentf90.out"))
    stop("seekparentf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(out, file.path(dir, "seekparentf90.out"))

  # Parse output files
  result <- list(output = out)

  # Check results (modified pedigree with match/no-match)
  ped_basename <- sub("\\.[^.]*$", "", basename(ped_file))
  check_file <- file.path(dir, paste0("Check_", basename(ped_file)))
  if (file.exists(check_file))
    result$check <- readLines(check_file)

  # Parent-offspring pair statistics
  pair_file <- file.path(dir, "Check_Parent_Pedigree.txt")
  if (file.exists(pair_file) && file.info(pair_file)$size > 0) {
    result$pair_checks <- readLines(pair_file)
  }

  # Sire search results
  sire_file <- file.path(dir, "Seek_Sire.txt")
  if (file.exists(sire_file) && file.info(sire_file)$size > 0) {
    result$seek_sire <- readLines(sire_file)
  }

  # Dam search results
  dam_file <- file.path(dir, "Seek_Dam.txt")
  if (file.exists(dam_file) && file.info(dam_file)$size > 0) {
    result$seek_dam <- readLines(dam_file)
  }

  # Assigned pedigree (corrected)
  assign_file <- file.path(dir, paste0("Assigned_", basename(ped_file)))
  if (file.exists(assign_file))
    result$assigned <- readLines(assign_file)

  return(result)
}
