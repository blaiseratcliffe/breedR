# Script which will be run as part of the installation used to automate install
# of BLUPF90+ binary. The script is run in a separate R environment containing
# the following variables: R_PACKAGE_NAME (the name of the package),
# R_PACKAGE_SOURCE (the path to the source directory of the package),
# R_PACKAGE_DIR (the path of the target installation directory of the package),
# R_ARCH (the arch-dependent part of the path, often empty), SHLIB_EXT (the
# extension of shared objects) and WINDOWS (TRUE on Windows, FALSE elsewhere).
# REF:
# http://cran.univ-paris1.fr/doc/manuals/r-release/R-exts.html#Package-subdirectories


source('../R/binaries.R')
source('../R/os.R')

# Set BREEDR_SKIP_INSTALL_BINARIES=true to install the package without
# downloading the BLUPF90+ backend (used in CI checks, where the binaries are
# not needed and may not run on the runner). Users can fetch them later with
# install_progsf90(). The download is also wrapped so a network failure does
# not abort installation.
if (identical(tolower(Sys.getenv("BREEDR_SKIP_INSTALL_BINARIES")), "true")) {
  message("BREEDR_SKIP_INSTALL_BINARIES is set; skipping BLUPF90+ download. ",
          "Run install_progsf90() to fetch the backend.")
} else {
  message("Downloading BLUPF90+ from:\n", breedr_progsf90_repo())
  tryCatch(
    install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin')),
    error = function(e)
      message("BLUPF90+ download failed (run install_progsf90() later): ",
              conditionMessage(e))
  )
}
