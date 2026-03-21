### RENUMF90 — Data preparation wrapper ###


#' Run RENUMF90 data preparation
#'
#' Runs the RENUMF90 program to renumber data, validate pedigrees, compute
#' inbreeding coefficients, and generate parameter files for BLUPF90+.
#'
#' RENUMF90 accepts alphanumeric IDs, handles unknown parent groups (UPGs),
#' computes inbreeding, and produces the \code{renf90.par} parameter file
#' that all downstream BLUPF90 programs consume.
#'
#' @param parfile path to a pre-written RENUMF90 parameter file. If provided,
#'   all other arguments except \code{dir} are ignored.
#' @param datafile path to the raw data file (space-separated, no header
#'   unless \code{skip_header} is used).
#' @param traits integer vector. Column positions of trait(s) in the data file.
#' @param residual_variance numeric matrix. Residual (co)variance matrix
#'   (n_traits x n_traits).
#' @param effects list of effect specifications. Each element is a list with:
#'   \describe{
#'     \item{pos}{integer vector of column positions (one per trait, 0 = missing)}
#'     \item{type}{character: "cross" or "cov"}
#'     \item{form}{character: "alpha" or "numer" (for cross-classified only)}
#'     \item{random}{character: "animal", "diagonal", or NULL for fixed}
#'     \item{optional}{character vector: "pe", "mat", "mpe" (animal effects only)}
#'     \item{file}{character: pedigree filename (animal/sire effects only)}
#'     \item{file_pos}{integer vector: animal, sire, dam, altdam, yob positions
#'       in pedigree file (default c(1,2,3,0,0))}
#'     \item{covariances}{numeric matrix: (co)variances for this random effect}
#'     \item{covariances_pe}{numeric matrix: PE (co)variances (if "pe" in optional)}
#'     \item{covariances_mpe}{numeric matrix: MPE (co)variances (if "mpe" in optional)}
#'     \item{nested}{list with pos and form for nested covariates}
#'   }
#' @param fields_passed integer vector. Columns passed through without
#'   renumbering (or NULL).
#' @param weights integer vector. Column position(s) for weights (or NULL).
#' @param snp_file character. Path to SNP genotype file (optional).
#' @param ped_depth integer. Pedigree search depth (0 = all, default 3).
#' @param inbreeding character. Inbreeding method: "pedigree" (default),
#'   "no-inbreeding", or "file <filename>".
#' @param upg_type character. Unknown parent group type: NULL (none), "yob",
#'   "in_pedigrees", "group", "group_unisex".
#' @param gen_int numeric vector of length 3: min, avg, max generation interval.
#' @param skip_header integer. Number of header lines to skip in data file.
#' @param options character vector. Additional OPTION lines for RENUMF90.
#' @param dir character. Working directory for input/output files.
#' @return A list with paths to output files and parsed statistics.
#' @export
renumf90 <- function(parfile = NULL,
                     datafile = NULL,
                     traits = NULL,
                     residual_variance = NULL,
                     effects = NULL,
                     fields_passed = NULL,
                     weights = NULL,
                     snp_file = NULL,
                     ped_depth = NULL,
                     inbreeding = NULL,
                     upg_type = NULL,
                     gen_int = NULL,
                     skip_header = NULL,
                     options = NULL,
                     dir = tempdir()) {

  bin_path <- breedR.getOption('breedR.bin')
  if (!check_renumf90(bin_path, quiet = TRUE))
    stop("renumf90 binary not installed. See ?install_renumf90", call. = FALSE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(parfile)) {
    # User provides pre-written parameter file
    if (!file.exists(parfile))
      stop("Parameter file not found: ", parfile, call. = FALSE)
    file.copy(parfile, file.path(dir, basename(parfile)), overwrite = TRUE)
    par_name <- basename(parfile)
  } else {
    # Build parameter file from R arguments
    if (is.null(datafile) || is.null(traits) || is.null(effects))
      stop("'datafile', 'traits', and 'effects' are required when ",
           "'parfile' is not provided.", call. = FALSE)

    par_name <- "renumf90.par"
    par_content <- build_renumf90_parfile(
      datafile = datafile, traits = traits,
      residual_variance = residual_variance, effects = effects,
      fields_passed = fields_passed, weights = weights,
      snp_file = snp_file, ped_depth = ped_depth,
      inbreeding = inbreeding, upg_type = upg_type,
      gen_int = gen_int, skip_header = skip_header,
      options = options
    )
    writeLines(par_content, file.path(dir, par_name))
  }

  # Run RENUMF90
  renum_out <- run_renumf90(dir, bin_path, par_name)
  writeLines(renum_out, file.path(dir, "renumf90.out"))

  # Parse outputs
  result <- parse_renumf90(dir)
  result$output <- renum_out
  result$dir <- dir

  return(result)
}


#' Build RENUMF90 parameter file content
#'
#' @inheritParams renumf90
#' @return Character vector of parameter file lines.
build_renumf90_parfile <- function(datafile, traits, residual_variance,
                                    effects, fields_passed, weights,
                                    snp_file, ped_depth, inbreeding,
                                    upg_type, gen_int, skip_header,
                                    options) {
  lines <- character(0)

  # DATAFILE
  lines <- c(lines, "DATAFILE", datafile)

  # SKIP_HEADER (optional)
  if (!is.null(skip_header))
    lines <- c(lines, "SKIP_HEADER", as.character(skip_header))

  # TRAITS
  lines <- c(lines, "TRAITS", paste(traits, collapse = " "))

  # FIELDS_PASSED TO OUTPUT
  lines <- c(lines, "FIELDS_PASSED TO OUTPUT")
  if (!is.null(fields_passed))
    lines <- c(lines, paste(fields_passed, collapse = " "))
  else
    lines <- c(lines, "")

  # WEIGHT(S)
  lines <- c(lines, "WEIGHT(S)")
  if (!is.null(weights))
    lines <- c(lines, paste(weights, collapse = " "))
  else
    lines <- c(lines, "")

  # RESIDUAL_VARIANCE
  lines <- c(lines, "RESIDUAL_VARIANCE")
  if (!is.null(residual_variance)) {
    rmat <- as.matrix(residual_variance)
    for (i in seq_len(nrow(rmat)))
      lines <- c(lines, paste(rmat[i, ], collapse = " "))
  } else {
    # Default: identity matrix scaled by 1
    n <- length(traits)
    lines <- c(lines, paste(rep(1, n), collapse = " "))
  }

  # EFFECTS
  for (eff in effects) {
    # Effect line: pos1 pos2 ... type [form]
    eff_line <- paste(c(eff$pos, eff$type), collapse = " ")
    if (eff$type == "cross" && !is.null(eff$form))
      eff_line <- paste(eff_line, eff$form)
    lines <- c(lines, "EFFECT", eff_line)

    # NESTED (optional, for covariates)
    if (!is.null(eff$nested)) {
      nest_line <- paste(c(eff$nested$pos, eff$nested$form), collapse = " ")
      lines <- c(lines, "NESTED", nest_line)
    }

    # RANDOM (if this effect is random)
    if (!is.null(eff$random)) {
      lines <- c(lines, "RANDOM", eff$random)

      # OPTIONAL (pe, mat, mpe — animal effects only)
      if (!is.null(eff$optional))
        lines <- c(lines, "OPTIONAL", paste(eff$optional, collapse = " "))

      # FILE (pedigree file for animal/sire)
      if (!is.null(eff$file))
        lines <- c(lines, "FILE", eff$file)

      # FILE_POS (pedigree column positions)
      if (!is.null(eff$file_pos))
        lines <- c(lines, "FILE_POS",
                    paste(eff$file_pos, collapse = " "))

      # SNP_FILE
      if (!is.null(snp_file) && eff$random == "animal")
        lines <- c(lines, "SNP_FILE", snp_file)

      # PED_DEPTH
      if (!is.null(ped_depth) && eff$random == "animal")
        lines <- c(lines, "PED_DEPTH", as.character(ped_depth))

      # GEN_INT
      if (!is.null(gen_int) && eff$random == "animal")
        lines <- c(lines, "GEN_INT", paste(gen_int, collapse = " "))

      # UPG_TYPE
      if (!is.null(upg_type) && eff$random == "animal")
        lines <- c(lines, "UPG_TYPE", upg_type)

      # INBREEDING
      if (!is.null(inbreeding) && eff$random == "animal")
        lines <- c(lines, "INBREEDING", inbreeding)

      # (CO)VARIANCES
      if (!is.null(eff$covariances)) {
        lines <- c(lines, "(CO)VARIANCES")
        cmat <- as.matrix(eff$covariances)
        for (i in seq_len(nrow(cmat)))
          lines <- c(lines, paste(cmat[i, ], collapse = " "))
      }

      # (CO)VARIANCES_PE
      if (!is.null(eff$covariances_pe)) {
        lines <- c(lines, "(CO)VARIANCES_PE")
        cmat <- as.matrix(eff$covariances_pe)
        for (i in seq_len(nrow(cmat)))
          lines <- c(lines, paste(cmat[i, ], collapse = " "))
      }

      # (CO)VARIANCES_MPE
      if (!is.null(eff$covariances_mpe)) {
        lines <- c(lines, "(CO)VARIANCES_MPE")
        cmat <- as.matrix(eff$covariances_mpe)
        for (i in seq_len(nrow(cmat)))
          lines <- c(lines, paste(cmat[i, ], collapse = " "))
      }
    }
  }

  # OPTION lines
  if (!is.null(options)) {
    for (opt in options)
      lines <- c(lines, paste("OPTION", opt))
  }

  return(lines)
}


## Execute renumf90 binary
run_renumf90 <- function(dir, bin_path, par_name = "renumf90.par") {

  renum_name <- renumf90_file(breedR.os.type())
  renum_src <- file.path(bin_path, renum_name)

  if (!file.exists(renum_src))
    stop("renumf90 binary not found at: ", renum_src,
         "\nInstall with install_renumf90().", call. = FALSE)

  # Copy binary to working directory (DLL isolation)
  renum_bin <- file.path(dir, renum_name)
  file.copy(renum_src, renum_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(renum_bin)
  })

  out <- system2(file.path(".", renum_name), input = par_name,
                 stdout = TRUE, stderr = TRUE)

  if (!is.null(attr(out, 'status'))) {
    writeLines(out, file.path(dir, "renumf90.out"))
    stop("renumf90 failed with exit code ",
         attr(out, 'status'), ".\nOutput:\n",
         paste(utils::tail(out, 30), collapse = "\n"),
         call. = FALSE)
  }

  return(out)
}


#' Parse RENUMF90 output files
#'
#' Reads the files produced by RENUMF90 and returns them as a structured list.
#'
#' @param dir directory containing RENUMF90 output files.
#' @return A list with components: par_file, data_file, ped_files, tables,
#'   inbreeding, xref_files.
parse_renumf90 <- function(dir) {

  result <- list()

  # Parameter file for BLUPF90
  par_file <- file.path(dir, "renf90.par")
  if (file.exists(par_file)) {
    result$par_file <- par_file
    result$par_content <- readLines(par_file)
  }

  # Renumbered data
  data_file <- file.path(dir, "renf90.dat")
  if (file.exists(data_file)) {
    result$data_file <- data_file
  }

  # Pedigree files (renaddXX.ped where XX is the effect number)
  ped_files <- list.files(dir, pattern = "^renadd.*\\.ped$", full.names = TRUE)
  if (length(ped_files) > 0) {
    result$ped_files <- ped_files
    # Parse the first (usually only) pedigree file
    result$pedigree <- utils::read.table(ped_files[1])
    ped_names <- c("animal", "parent1", "parent2", "inb_upg_code",
                    "yob", "known_parents", "n_records",
                    "n_progeny_p1", "n_progeny_p2", "original_id")
    if (ncol(result$pedigree) >= 10)
      names(result$pedigree) <- ped_names[seq_len(ncol(result$pedigree))]
  }

  # Renumbering tables
  tables_file <- file.path(dir, "renf90.tables")
  if (file.exists(tables_file)) {
    result$tables_file <- tables_file
    result$tables_content <- readLines(tables_file)
  }

  # Field descriptions
  fields_file <- file.path(dir, "renf90.fields")
  if (file.exists(fields_file)) {
    result$fields_file <- fields_file
  }

  # Inbreeding coefficients
  inb_file <- file.path(dir, "renf90.inb")
  if (file.exists(inb_file) && file.info(inb_file)$size > 0) {
    result$inbreeding <- utils::read.table(inb_file)
    names(result$inbreeding) <- c("original_id", "inbreeding")[
      seq_len(ncol(result$inbreeding))]
  }

  # XrefID files (for genomic analyses)
  xref_files <- list.files(dir, pattern = "_XrefID$", full.names = TRUE)
  if (length(xref_files) > 0) {
    result$xref_files <- xref_files
  }

  return(result)
}


#' Fit a model using RENUMF90 output
#'
#' Takes the output of \code{\link{renumf90}} and runs BLUPF90+ on the
#' renumbered data and parameter file. This bridges RENUMF90 data preparation
#' with BLUPF90+ model fitting.
#'
#' This function is particularly useful for multi-trait models with different
#' effects per trait, which the \code{\link{remlf90}} formula interface does
#' not support. Use \code{\link{renumf90}} to set up the model with per-trait
#' effect positions, then fit with this function.
#'
#' @param renum output from \code{\link{renumf90}}.
#' @param method either 'ai' or 'em'.
#' @param progsf90.options character vector. Additional OPTIONS for BLUPF90+.
#'   Use \code{\link{var_functions}}, \code{\link{h2_formula}}, or
#'   \code{\link{rg_formula}} to compute genetic parameters with SEs.
#' @param debug logical.
#' @return A list with model results (solutions, variance components, etc.)
#'   and a \code{renum} component containing the RENUMF90 outputs.
#'
#' @examples
#' \dontrun{
#' ## Multi-trait model: height and diameter with different fixed effects
#' ## Data file columns: 1=id, 2=herd, 3=age, 4=site, 5=height, 6=diameter
#'
#' renum <- renumf90(
#'   datafile = "data.txt",
#'   traits = c(5, 6),
#'   residual_variance = matrix(c(10, 3, 3, 5), 2, 2),
#'   effects = list(
#'     # Herd: both traits (column 2)
#'     list(pos = c(2, 2), type = "cross", form = "alpha"),
#'     # Age: height only (column 3, 0 = absent for diameter)
#'     list(pos = c(3, 0), type = "cov"),
#'     # Site: diameter only (0 = absent for height, column 4)
#'     list(pos = c(0, 4), type = "cross", form = "alpha"),
#'     # Genetic effect: both traits
#'     list(pos = c(1, 1), type = "cross", form = "alpha",
#'          random = "animal", file = "pedigree.txt",
#'          covariances = matrix(c(5, 2, 2, 3), 2, 2))
#'   )
#' )
#'
#' ## Fit with heritabilities and genetic correlation
#' res <- remlf90_from_renum(renum, method = 'ai',
#'   progsf90.options = c(
#'     h2_formula(genetic_effect = 4, trait = 1),
#'     h2_formula(genetic_effect = 4, trait = 2),
#'     rg_formula(trait1 = 1, trait2 = 2, effect = 4)
#'   )
#' )
#' }
#' @export
remlf90_from_renum <- function(renum,
                                method = c('ai', 'em'),
                                progsf90.options = NULL,
                                debug = FALSE) {

  method <- match.arg(method)

  if (is.null(renum$par_file) || !file.exists(renum$par_file))
    stop("RENUMF90 parameter file not found. Run renumf90() first.",
         call. = FALSE)

  dir <- renum$dir
  bin_path <- breedR.getOption('breedR.bin')

  if (!check_progsf90(bin_path, quiet = TRUE))
    stop("BLUPF90+ binary not installed. See ?install_progsf90",
         call. = FALSE)

  # Read the RENUMF90-generated parameter file
  par_lines <- readLines(renum$par_file)

  # Add method-specific options
  method_opts <- 'method VCE'
  if (method == 'em') method_opts <- c(method_opts, 'EM-REML')
  par_lines <- c(par_lines,
                  paste("OPTION", c("sol se", method_opts)),
                  if (!is.null(progsf90.options))
                    paste("OPTION", progsf90.options))

  # Write the updated parameter file
  writeLines(par_lines, file.path(dir, "renf90.par"))

  # Execute BLUPF90+
  blupf90_name <- progsf90_files(breedR.os.type())
  blupf90_src <- file.path(bin_path, blupf90_name)

  # Copy binary to working directory (DLL isolation)
  blupf90_bin <- file.path(dir, blupf90_name)
  file.copy(blupf90_src, blupf90_bin, overwrite = TRUE)

  cdir <- setwd(dir)
  on.exit({
    setwd(cdir)
    unlink(blupf90_bin)
  })

  if (!debug) {
    reml.out <- system2(file.path(".", blupf90_name),
                        input = 'renf90.par',
                        stdout = TRUE, stderr = TRUE)
  } else {
    cat("Parameter file:\n")
    cat(par_lines, sep = "\n")
    return(invisible(NULL))
  }

  if (!is.null(attr(reml.out, 'status'))) {
    stop("BLUPF90+ failed with exit code ",
         attr(reml.out, 'status'), ".\nOutput:\n",
         paste(utils::tail(reml.out, 20), collapse = "\n"),
         call. = FALSE)
  }

  writeLines(reml.out, file.path(dir, "blupf90.out"))

  # Parse solutions
  sol_file <- file.path(dir, "solutions")
  result <- list()

  if (file.exists(sol_file)) {
    sol <- utils::read.table(sol_file, header = TRUE, row.names = NULL)
    names(sol) <- c("trait", "effect", "level", "solution",
                     if (ncol(sol) >= 5) "se")[seq_len(ncol(sol))]
    result$solutions <- sol
  }

  result$output <- reml.out
  result$method <- method
  result$renum <- renum

  return(result)
}
