## Internal utility functions
## Not exported

# lme4-style formulas
# 
# Transform the separated fixed and random formulas
# into the single formula with lme4 syntaxis
lme4_fml <- function(fix, rnd, rm_int = TRUE) {
  rnd.terms <- attr(stats::terms(rnd), 'term.labels')
  rnd.terms.lme4 <- paste('(1|', rnd.terms, ')', sep ='')
  int <- ifelse(rm_int, '-1', '')
  rnd.upd <- paste('~ .', int, paste('+', rnd.terms.lme4, collapse = ' '))
  fml.lme4 <- update(fix, rnd.upd)
  return(fml.lme4)
}


breedR.is.element <- function(name, alist) {
  ## return TRUE if element with name NAME is a member of LIST and
  ## the value is non null and not NA.
  if (any(names(alist) == name)) {
    idx = which(names(alist) == name)
    if (!is.null(alist[[idx]])) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  } else {
    return (FALSE)
  }
}

breedR.get.element <-  function(name, alist) {
  if (breedR.is.element(name, alist)) {
    return (alist[[which(names(alist) == name)]])
  } else {
    return (NULL)
  }
}


# Geometric mean
gmean <- function(x) {
  logx <- log(x)
  finite.logx <- is.finite(logx)
  if( !all(finite.logx) ) {
    warning('Removing zeroes for geometric mean')
    logx <- logx[finite.logx]
  }
  
  return(exp(mean(logx)))
}


# BreedR binaries
# 
# Return the default path to breedR binaries.
# All BLUPF90+ binaries are stored flat in the 'bin' directory.
breedR.bin.builtin  <- function()
{
  if (breedR.os.type() != 'else') {
    return(file.path(system.file(package = 'breedR'), 'bin'))
  } else {
    stop("Unknown platform")
  }
}


#' Determine the user's home directory
#' 
#' Relies on \code{Sys.getenv('HOME')}, or under windows, on
#' \code{Sys.getenv("USERPROFILE"))} changing backslashes to slashes.
`breedR.get.HOME` = function()
{
  return (as.character(ifelse(breedR.os("windows"),
                              gsub("\\\\", "/", Sys.getenv("USERPROFILE")),
                              Sys.getenv("HOME"))))
}

#' Determine the user name
`breedR.get.USER` = function()
{
  u = ""
  for (U in c("USER", "USERNAME", "LOGNAME")) {
    u = Sys.getenv(U)
    if (u != "")
      break;
  }
  if (u == "")
    u = "UnknownUserName"
  
  return (as.character(u))
}

# Convert a incidence matrix specified in 8+8 columns format
# the first 8 are coefficients, and the last 8 are columns
# into a sparse matrix format
matrix.short16 <- function(M) {
  coef = M[, 1:8]
  neig = M[, 8+1:8]
  
  n <- nrow(coef)
  p <- max(neig, na.rm = TRUE)
  
  i <- rep(1:n, 8)
  j <- as.vector(neig)
  x <- as.vector(coef)
  
  rm.idx <- which(x==0)
  stopifnot( all(x[rm.idx] == 0) )
  
  Z <- Matrix::spMatrix(nrow = n, ncol = p,
                        i = i[-rm.idx],
                        j = j[-rm.idx],
                        x = x[-rm.idx])
  return(Z)
}

#' 'Splat' arguments to a function
#' 
#' Wraps a function in do.call, so instead of taking multiple arguments, it
#' takes a single named list which will be interpreted as its arguments.
#' 
#' @param flat function to splat
#' 
#' This is useful when you want to pass a function a row of data frame or array,
#' and don't want to manually pull it apart in your function.
#' 
#' Borrowed from \code{\link[plyr]{splat}}
#' 
#' @return a function
#' 
#' @examples 
#'   args <- replicate(3, runif(5), simplify = FALSE)
#'   identical(breedR:::splat(rbind)(args), do.call(rbind, args))
splat <- function (flat) {
  function(args, ...) {
    do.call(flat, c(args, list(...)))
  }
}

# given a list of data.frames, extract a given column
# from each of the data.frames into a matrix.
# Optionally drop into a vector if dimension = 1.
# Used in ranef.remlf90 and fixef.remlf90 to extract
# trait-wise predictions of effects
ldf2matrix <- function(x, vname, drop = TRUE) {
  ## All dataframes (ntraits) are of the same size (nlevels x (value, s.e.))
  ## ensure a matrix, even if nlevels = 1
  ans <- do.call(cbind, lapply(x, `[[`, vname))
  rownames(ans) <- rownames(x[[1]])
  if (drop) ans <- drop(ans)
  return(ans)
}


# Extract values and standard errors from lists
# of effects estimates (or predictions)
# x is a list of effects, where each element is a trait-wise
# list of data.frames with columns 'value' and 's.e.'
get_estimates <- function(x) {
  values <- lapply(x, ldf2matrix, 'value')
  se <- lapply(x, ldf2matrix, 's.e.')
  ans <- mapply(function(gvl, gse) structure(gvl, se = gse),
                values, se, SIMPLIFY = FALSE)
  return(ans)
}

# combine sub-effect names and trait names
# trait names within sub-effect names
# bl1 bl2 bl3 + y1 y2 = bl1.y1 bl1.y2 bl2.y1 bl2.y2 bl3.y1 bl3.y2
# names_effect(paste0("bl", 1:3), paste0("y", 1:2))
# "bl1.y1" "bl1.y2" "bl2.y1" "bl2.y2" "bl3.y1" "bl3.y2"
# names_effect(paste0("bl", 1:3), NULL)
# "bl1" "bl2" "bl3"
# names_effect(NULL, paste0("y", 1:2))
# "y1" "y2"
# names_effect(NULL, NULL)
# NULL
names_effect <- function(inner = NULL, outer = NULL) {

  ans <- inner
  if (length(outer) > 1) {
    if (!is.null(inner)) {
      ans <- apply(
        expand.grid(
          outer, inner,
          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
        )[2:1],
        1, paste, collapse = "."
      )
    } else {
      ans <- outer
    }
  }
  return(ans)
}


## component names
vcnames <- function(efname, efdim, trnames) {
  
  dim_subtrait <- efdim/ifelse(is.null(trnames), 1, length(trnames))
  if (dim_subtrait > 1)
    efname <- paste(efname, seq_len(dim_subtrait), sep = "_")
  diag_names <- names_effect(efname, trnames)
  
  ## matrix components include variances and covariances
  ans <- outer(diag_names, diag_names, paste, sep = "_")
  diag(ans) <- diag_names
  return(t(ans)[lower.tri(ans, diag = TRUE)])
}


## transform a list of 2 cov-matrices (est; SE) into a 2-col data frame
## with properly named estimates. Use nm to name the original effect
lmat2df <- function(x, nm) {
  data.frame(
    lapply(x, `[`, lower.tri(x[[1]], diag = TRUE)),
    row.names =  vcnames(nm, nrow(x[[1]]), trnames = rownames(x[[1]])),
    check.names = FALSE
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Working directories     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# A fresh working directory under tempdir().
#
# Every run that drives a PROGSF90 binary needs one of its own: the programs
# write a fixed set of names (parameters, data, solutions), so two runs sharing
# a directory means the second silently consumes the first one's files. Used by
# remlf90(), gibbsf90() and breedR.qget() alike.
#
# Warnings from dir.create() are deliberately not suppressed: the reason it
# failed -- permissions, a full disk -- is exactly what the caller's error
# message cannot say on its own.
breedR_workdir <- function(prefix = "breedR_") {
  d <- tempfile(prefix, tmpdir = tempdir())
  dir.create(d, recursive = TRUE)
  if (!dir.exists(d))
    stop("Could not create the working directory: ", d, call. = FALSE)
  d
}


#' Remove the working directory of a fit
#'
#' Each fit or sampling run works in its own directory under \code{tempdir()},
#' holding its parameter file, data, structure files and solutions. They are
#' kept so that follow-up steps such as \code{\link{postgsf90}} can read them,
#' and are removed when the session ends. A session that fits many models keeps
#' one directory per fit, which is worth reclaiming explicitly when the models
#' are large or numerous.
#'
#' Only directories under \code{tempdir()} are removed. \code{\link{renumf90}}
#' and friends accept a user-chosen \code{dir}, and deleting one of those on the
#' strength of a stored path is not a mistake worth risking, so it is refused.
#'
#' @param x a fitted model from \code{\link{remlf90}}, a result from
#'   \code{\link{gibbsf90}}, \code{\link{postgsf90}} or \code{\link{renumf90}},
#'   or a character path to the directory itself.
#' @return \code{TRUE} if a directory was removed, \code{FALSE} if there was
#'   nothing to remove (no recorded directory, or already gone). Invisibly.
#' @seealso \code{\link{remlf90}} for \code{res$reml$dir}.
#' @examples
#' \dontrun{
#'   res <- remlf90(phe_X ~ gg, data = globulus)
#'   clean_workdir(res)
#' }
#' @export
clean_workdir <- function(x) {

  dir <- if (is.character(x)) x
         else if (!is.null(x$reml$dir)) x$reml$dir
         else x$dir

  if (is.null(dir) || !nzchar(dir)) return(invisible(FALSE))
  if (!dir.exists(dir)) return(invisible(FALSE))

  ## Compare resolved paths: the stored one and tempdir() can disagree on
  ## separators and short names on Windows while naming the same place.
  full <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)

  ## Strictly *under* tempdir(). Equal to it is the session's own scratch
  ## space, whose removal breaks everything downstream of here, so a stored
  ## path that resolves to bare tempdir() is refused rather than obeyed.
  if (identical(full, root) || !startsWith(full, paste0(root, "/")))
    stop("Refusing to remove ", full, ", which is not a directory under ",
         "tempdir().\n",
         "  Remove it yourself if that is really what you want.",
         call. = FALSE)

  unlink(full, recursive = TRUE)
  invisible(!dir.exists(full))
}
