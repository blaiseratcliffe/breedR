### Functions to download and install BLUPF90 binary dependencies ###

## UGA platform directory names for download URLs
uga_platform_dir <- function(platform = breedR.os.type()) {
  switch(platform,
         windows = "Windows",
         linux = "Linux",
         mac = "Mac_OSX",
         stop("Unsupported platform: ", platform))
}

#' Default repository for BLUPF90 binaries
#'
#' Returns the UGA download URL. Can be overridden with the
#' PROGSF90_URL environment variable.
breedr_progsf90_repo <- function() {
  if (!nchar(url <- Sys.getenv("PROGSF90_URL"))) {
    url <- "https://nce.ads.uga.edu/html/projects/programs"
  }
  return(url)
}


## Return the expected filename of the BLUPF90+ binary by platform.
## BLUPF90+ replaces the standalone remlf90 and airemlf90 programs.
progsf90_files <- function(os = breedR.os.type(),
                           compressed = FALSE) {
  if (compressed) return(character(0))  # no longer using archives
  ans <- "blupf90+"
  if (os == 'windows') ans <- paste0(ans, ".exe")
  return(ans)
}


#' Checks installation of BLUPF90+ binary
#'
#' Checks whether the BLUPF90+ binary is installed in the right directory.
#' If not, allows calling the installer.
#'
#' @param path directory to check for the presence of binaries. Default is
#'   defined in the package options, and it depends on the platform.
#' @param platform Either "linux", "windows" or "mac". Default is
#'   current platform.
#' @param quiet if TRUE, it won't ask whether to install missing binaries.
#' @export
check_progsf90 <- function(path = breedR.getOption('breedR.bin'),
                           platform = breedR.os.type(),
                           quiet = !interactive() ) {

  bin.list <- progsf90_files(platform, compressed = FALSE)

  check <- FALSE
  if (file.exists(path)) {
    check <- all(bin.list %in% list.files(path))
  }

  if (!check && !quiet) {
    message('Binary dependencies missing.',
            '\nWould you like to install them?\t')
    if (utils::menu(c("Yes", "No")) == 1) {
      install_progsf90(dest = path, platform = platform)
      check <- check_progsf90(path, platform, quiet)
    }
  }

  return(invisible(check))
}

#' Install BLUPF90+ binary
#'
#' Downloads the BLUPF90+ unified binary from the UGA BLUPF90 distribution.
#'
#' @param url base URL for the BLUPF90 program repository.
#' @param dest destination directory for the binary. Default is 'bin' under
#'   the current installation dir.
#' @param platform what version of the binary to install. Default is
#'   current platform.
#' @param arch Either "32" or "64". Coerced to string if necessary.
#' @param quiet logical. Whether not to display messages.
#' @export
install_progsf90 <- function(
  url      = breedr_progsf90_repo(),
  dest     = system.file('bin', package = 'breedR'),
  platform = breedR.os.type(),
  arch     = breedR.os.32or64bit(),
  quiet    = !interactive()
) {
  execs <- progsf90_files(platform, compressed = FALSE)
  f.url <- file.path(url, uga_platform_dir(platform), paste0(arch, "bit"))

  res <- vapply(execs, function(f) {
    retrieve_bin_direct(f, url = f.url, dest = dest)
  }, logical(1))

  return(res)
}


### Auxiliary program binary management ###

## Return expected filenames for auxiliary programs by platform.
## Includes genomic programs (preGSf90, postGSf90) and data
## preparation (renumf90).
genomic_program_files <- function(os = breedR.os.type()) {
  progs <- c("preGSf90", "postGSf90", "predf90",
             "validationf90", "predictf90", "seekparentf90",
             "gibbsf90+", "postgibbsf90", "qcf90")
  if (os == 'windows') progs <- paste0(progs, ".exe")
  return(progs)
}

## Return expected filename for renumf90 by platform
renumf90_file <- function(os = breedR.os.type()) {
  f <- "renumf90"
  if (os == 'windows') f <- paste0(f, ".exe")
  return(f)
}

#' Check installation of renumf90 binary
#'
#' @param path directory to check for the binary.
#' @param platform Either "linux", "windows" or "mac".
#' @param quiet if TRUE, won't prompt to install.
#' @return Logical (invisible).
#' @export
check_renumf90 <- function(
  path = breedR.getOption('breedR.bin'),
  platform = breedR.os.type(),
  quiet = !interactive()
) {
  bin <- renumf90_file(platform)
  check <- file.exists(path) && bin %in% list.files(path)

  if (!check && !quiet) {
    message('renumf90 binary missing.\nWould you like to install it?\t')
    if (utils::menu(c("Yes", "No")) == 1) {
      install_renumf90(dest = path, platform = platform)
      check <- check_renumf90(path, platform, quiet = TRUE)
    }
  }
  return(invisible(check))
}

#' Install renumf90 binary
#'
#' Downloads renumf90 from the UGA BLUPF90 distribution.
#'
#' @param url base URL for the BLUPF90 program repository.
#' @param dest destination directory.
#' @param platform Either "linux", "windows", or "mac".
#' @param arch Either "32" or "64".
#' @param quiet logical.
#' @return Logical indicating success.
#' @export
install_renumf90 <- function(
  url      = breedr_progsf90_repo(),
  dest     = breedR.getOption('breedR.bin'),
  platform = breedR.os.type(),
  arch     = breedR.os.32or64bit(),
  quiet    = !interactive()
) {
  f <- renumf90_file(platform)
  f.url <- file.path(url, uga_platform_dir(platform), paste0(arch, "bit"))
  retrieve_bin_direct(f, url = f.url, dest = dest)
}

#' Check installation of genomic program binaries
#'
#' Checks whether preGSf90 and postGSf90 are installed. If not, allows
#' calling the installer. These programs are optional and only needed for
#' genomic selection (ssGBLUP) and genome-wide association (ssGWAS).
#'
#' @param path directory to check for the presence of binaries. Default is
#'   defined in the package options.
#' @param platform Either "linux", "windows" or "mac". Default is
#'   current platform.
#' @param quiet if TRUE, it won't ask whether to install missing binaries.
#' @return Logical (invisible) indicating whether binaries are present.
#' @export
check_genomic_programs <- function(
  path = breedR.getOption('breedR.bin'),
  platform = breedR.os.type(),
  quiet = !interactive()
) {
  bin.list <- genomic_program_files(platform)

  check <- FALSE
  if (file.exists(path)) {
    check <- all(bin.list %in% list.files(path))
  }

  if (!check && !quiet) {
    message('Genomic program binaries (preGSf90, postGSf90) missing.',
            '\nWould you like to install them?\t')
    if (utils::menu(c("Yes", "No")) == 1) {
      install_genomic_programs(dest = path, platform = platform)
      check <- check_genomic_programs(path, platform, quiet = TRUE)
    }
  }

  return(invisible(check))
}

#' Install genomic program binaries
#'
#' Downloads preGSf90 and postGSf90 from the UGA BLUPF90 distribution.
#'
#' @param url base URL for the BLUPF90 program repository.
#' @param dest destination directory for the binaries.
#' @param platform Either "linux", "windows", or "mac". Default is
#'   current platform.
#' @param arch Either "32" or "64".
#' @param quiet logical. Whether not to display messages.
#' @return Named logical vector indicating success for each program.
#' @export
install_genomic_programs <- function(
  url      = breedr_progsf90_repo(),
  dest     = breedR.getOption('breedR.bin'),
  platform = breedR.os.type(),
  arch     = breedR.os.32or64bit(),
  quiet    = !interactive()
) {
  execs <- genomic_program_files(platform)
  f.url <- file.path(url, uga_platform_dir(platform), paste0(arch, "bit"))

  res <- vapply(execs, function(f) {
    retrieve_bin_direct(f, url = f.url, dest = dest)
  }, logical(1))

  return(res)
}


## Download a raw executable (no decompression needed).
## Sets execute permissions on Unix platforms.
retrieve_bin_direct <- function(f, url, dest) {
  destf <- file.path(dest, f)
  if (!file.exists(dest))
    dir.create(dest, recursive = TRUE)

  out <- tryCatch(
    utils::download.file(
      url = file.path(url, f),
      destfile = destf,
      mode = 'wb',
      cacheOK = FALSE,
      quiet = TRUE
    ),
    error = identity
  )

  if (inherits(out, 'error')) {
    warning("Download of ", f, " failed: ", conditionMessage(out))
    unlink(destf)
    return(FALSE)
  }

  # Set execute permissions on Unix
  if (.Platform$OS.type == "unix") {
    Sys.chmod(destf, mode = "0755")
  }

  return(file.exists(destf) && file.info(destf)$size > 0)
}
