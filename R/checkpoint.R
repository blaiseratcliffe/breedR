## Checkpointing and resuming a REML fit.
##
## AIREMLF90 prints the full state of the optimiser -- the residual and genetic
## (co)variance matrices -- at every round. A streamed log is therefore a
## complete per-round checkpoint, and resuming a fit is nothing more than
## reading the last usable round back into the initial variances.
##
## The blocks are laid out as
##
##     In round <N>  convergence= ...
##     delta convergence= ...
##    new R
##     <residual (co)variance matrix>
##    new G
##     <(co)variance matrix of the first random group>
##    new G
##     ...
##
## and the `new G` blocks appear in the same order as the
## `Genetic variance(s) for effect K` blocks of the final section, which is the
## order of names(effects)[effect_type == 'random'].
##
## Note that the residual comes *first* within a round but *last* in the final
## estimates; the two representations are not positionally interchangeable.


#' Line numbers of the round headers of a REML log
#'
#' @param x character vector. Lines of a REML log.
#'
#' @return integer vector of line numbers, named with the round numbers.
#'   \code{integer(0)} if the log contains no completed round.
#' @keywords internal
reml_round_index <- function(x) {

  re <- "^[[:space:]]*In round[[:space:]]+([0-9]+)"
  idx <- grep(re, x)

  structure(idx,
            names = sub(paste0(re, ".*$"), "\\1", x[idx]))
}


#' Extract a (co)variance block, refusing blocks that run to the end of file
#'
#' A log that was truncated by a crash ends in the middle of a line rather than
#' on a line boundary. \code{parse.txtmat} does not always notice: a value cut
#' inside its exponent, say \code{0.57309E+09} truncated to \code{0.573}, parses
#' cleanly and silently rescales a variance component by nine orders of
#' magnitude. No amount of error trapping catches that, because nothing goes
#' wrong.
#'
#' The invariant that does catch it: a block which is genuinely finished is
#' followed by more text, so a block reaching the last line of the log is never
#' trustworthy. This also makes the parser correct on a log that is still being
#' written, where the final block is always incomplete.
#'
#' @param l numeric. Line where the block starts (just after the heading).
#' @param x character vector. Lines of the log.
#' @param text_lines integer vector. Non-numeric line numbers of \code{x},
#'   computed once by the caller. Must come from this same \code{x}.
#'
#' @return numeric matrix, or \code{NULL} if the block is unusable.
#' @keywords internal
parse_block_strict <- function(l, x, text_lines = which(!is_numericlog(x))) {

  tryCatch({
    blk <- extract_block(l, x, text_lines)

    ## Refuse a block that is not terminated by a following line.
    if (l + length(blk) - 1L >= length(x)) return(NULL)

    ## Wrapped matrices warn ("data length is not a multiple of split
    ## variable") on some truncations; the value is unusable either way.
    suppressWarnings(parse.txtmat(blk))
  },
  error = function(e) NULL)
}


#' (Co)variance matrices printed within one round of a REML log
#'
#' Three dialects are in circulation and all three must be read: BLUPF90+ under
#' AI-REML emits \code{new R} / \code{new G}, the same binary under EM-REML emits
#' \code{new r} / \code{new G}, and the legacy REMLF90 emits \code{new r} and a
#' bare \code{G}. Matching only one pair silently recovers nothing from the
#' others.
#'
#' @param x character vector. The whole log.
#' @param from,to integer. Bounds of the round segment, inclusive.
#' @param text_lines integer vector. Non-numeric line numbers of \code{x},
#'   computed once by the caller. Must come from this same \code{x}.
#'
#' @return list with elements \code{R} (matrix or \code{NULL}) and \code{G}
#'   (list of matrices, in parameter-file order).
#' @keywords internal
reml_round_covariances <- function(x, from, to,
                                   text_lines = which(!is_numericlog(x))) {

  ## Anchored on the whole log so that extract_block() sees the real
  ## surroundings of the block, and never the parameter echo near the top of
  ## the file -- parse.txtmat() drops the first field of every line, which is
  ## only safe where Fortran carriage control guarantees a blank first column.
  res_re <- "^[[:space:]]*new [Rr][[:space:]]*$"
  grp_re <- "^[[:space:]]*(new [Gg]|G)[[:space:]]*$"

  win <- seq.int(from, min(to, length(x)))

  r_at <- win[grepl(res_re, x[win])]
  g_at <- win[grepl(grp_re, x[win])]

  list(R = if (length(r_at))
             parse_block_strict(utils::tail(r_at, 1) + 1L, x, text_lines)
           else NULL,
       G = lapply(g_at + 1L, parse_block_strict, x = x,
                  text_lines = text_lines))
}


#' Last usable checkpoint in a REML log
#'
#' Walks backwards from the final round and returns the first one that is
#' complete, conformant with the model at hand and (optionally) usable as a set
#' of initial variances.
#'
#' The positive-definiteness test is applied to values the backend printed at
#' about five significant figures, so it can disagree with the backend in either
#' direction. It is still worth having: without it a resumed fit can be seeded
#' from a state the optimiser had already rejected, which
#' \code{validate_variance} would then refuse anyway.
#'
#' @param x character vector. Lines of a REML log.
#' @param n_groups integer. Number of random-effect groups expected.
#' @param dims integer vector of length \code{n_groups + 1}. Expected dimension
#'   of each group and, last, of the residual. \code{NULL} to skip the check.
#' @param names character vector of length \code{n_groups + 1}, ending in
#'   \code{'residuals'}. \code{NULL} for positional names.
#' @param spd logical. Require every recovered matrix to be positive definite.
#'
#' @return named list of matrices with an integer attribute \code{round}. When
#'   no round qualifies, an empty list with a character attribute
#'   \code{reason}; test for failure with \code{!length(.)}.
#' @keywords internal
last_reml_checkpoint <- function(x, n_groups, dims = NULL,
                                 names = NULL, spd = TRUE) {

  ## NULL cannot carry attributes, so failure is an empty list.
  failed <- function(reason) structure(list(), reason = reason)

  anchors <- reml_round_index(x)
  if (!length(anchors))
    return(failed("no completed round"))

  ## Classify every line once. extract_block() would otherwise rescan the whole
  ## log for each block, which is one scan per random group per round: on a
  ## 542-round, 16,686-line log a full walk back took 35 s, essentially all of
  ## it rescanning.
  ##
  ## This must stay inside this function, keyed to the `x` it was given. Hoist
  ## it any higher and it could be paired with a different vector -- notably
  ## reml_checkpoint_var.ini(), which chooses between a caller-supplied `lines`
  ## and a fresh read of a file that may still be growing. A stale index gives
  ## block bounds past the end of `x`, and parse_block_strict() would swallow
  ## the resulting error as NULL: a silently skipped round, not a failure.
  text_lines <- which(!is_numericlog(x))

  if (is.null(names))
    names <- c(if (n_groups) paste0("G", seq_len(n_groups)), "residuals")

  is_spd <- function(m) {
    isSymmetric(unname(m), check.attributes = FALSE) &&
      all(eigen(m, symmetric = TRUE, only.values = TRUE)$values > 0)
  }

  reason <- "no completed round"

  for (i in rev(seq_along(anchors))) {

    from <- anchors[[i]]
    to <- if (i < length(anchors)) anchors[[i + 1L]] - 1L else length(x)

    cv <- reml_round_covariances(x, from, to, text_lines)
    blocks <- c(cv$G, list(cv$R))

    ## A round the backend never finished writing, or that a crash cut short.
    if (is.null(cv$R) || any(vapply(cv$G, is.null, TRUE))) {
      if (i == length(anchors)) reason <- "incomplete"
      next
    }

    ## A different model: it printed a different number of groups. Reported
    ## separately from incompleteness, because the remedy is different.
    if (length(cv$G) != n_groups) {
      reason <- "dimension"
      next
    }

    ## Conformant with this model?
    if (!is.null(dims) &&
        !all(vapply(seq_along(blocks),
                    function(j) all(dim(blocks[[j]]) == dims[[j]]), TRUE))) {
      reason <- "dimension"
      next
    }

    ## Usable as initial variances?
    if (spd && !all(vapply(blocks, is_spd, TRUE))) {
      reason <- "not positive definite"
      next
    }

    return(structure(stats::setNames(blocks, names),
                     round = as.integer(base::names(anchors)[i])))
  }

  failed(reason)
}


#' Random-effect names and dimensions of a built model
#'
#' The key property this whole feature rests on: the vector returned here is the
#' one \code{parse_results} uses to label the variance components, so a
#' checkpoint written by one run is keyed identically when read by the next.
#'
#' @param effects list of effects, as built by \code{build.effects}.
#' @param ntraits integer.
#'
#' @return list with \code{names} and \code{dims}.
#' @keywords internal
reml_checkpoint_layout <- function(effects, ntraits) {

  rand <- which(vapply(effects, effect_type, '') == 'random')

  list(names = c(names(effects)[rand], 'residuals'),
       dims = c(vapply(effects[rand],
                       function(g) nrow(as.matrix(g$cov.ini)), 1L),
                ntraits),
       n_groups = length(rand))
}


#' Recover initial variances from a previous run's log
#'
#' @param file character. Path of the log, used for messages.
#' @param effects list of effects, as built by \code{build.effects}.
#' @param ntraits integer.
#' @param lines character vector. Lines of the log, when already read.
#' @param spd logical. Passed to \code{last_reml_checkpoint}.
#'
#' @return list with \code{var} (named list of matrices) and \code{round}.
#' @keywords internal
reml_checkpoint_var.ini <- function(file, effects, ntraits,
                                    lines = NULL, spd = TRUE) {

  x <- if (is.null(lines)) readLines(file, warn = FALSE) else lines
  lay <- reml_checkpoint_layout(effects, ntraits)

  ck <- last_reml_checkpoint(x, n_groups = lay$n_groups, dims = lay$dims,
                             names = lay$names, spd = spd)

  if (!length(ck)) {
    reason <- attr(ck, 'reason')
    detail <-
      switch(reason,
             dimension = paste0(
               "its rounds do not match this model, which has ",
               lay$n_groups, " random-effect group(s) (",
               paste(utils::head(lay$names, -1), collapse = ", "),
               "). The log was probably written by a different model"),
             `not positive definite` = paste(
               "none of its rounds holds a positive-definite set of",
               "(co)variance matrices"),
             incomplete = "it holds no complete round",
             paste0("it holds ", reason))
    stop("Cannot resume from '", file, "': ", detail, ".", call. = FALSE)
  }

  ## The constructor's validation is bypassed by assigning cov.ini directly,
  ## so run it here instead.
  for (nm in utils::head(lay$names, -1))
    validate_variance(ck[[nm]],
                      dimension = dim(as.matrix(effects[[nm]]$cov.ini)),
                      what = paste0("recovered variance for '", nm, "'"),
                      where = paste("progress file", file))

  list(var = ck, round = attr(ck, 'round'))
}


#' Read the variance components of a breedR REML log
#'
#' Recovers the (co)variance matrices of the last usable round of a REML log,
#' including one that is still being written. This makes the current state of a
#' long-running fit visible without interrupting it, and lets the components of
#' a crashed or non-converged run be reused.
#'
#' Because \code{remlf90} blocks while the backend runs, checking on a fit in
#' progress means calling this from a \emph{second} R session.
#'
#' The last round of the log is deliberately never used: a log truncated by a
#' crash ends mid-line, and a value cut inside its exponent parses cleanly while
#' being wrong by orders of magnitude. The last round that is followed by
#' further output is the last one that can be trusted.
#'
#' @param file character. Path of a REML log, e.g. one written by the
#'   \code{progress_file} argument of \code{\link{remlf90}}.
#' @param model a fitted \code{remlf90} object. When given, the components are
#'   named after its random effects and checked against their dimensions.
#'   Otherwise names are positional.
#' @param spd logical. Only return a round whose matrices are all positive
#'   definite. \code{FALSE} returns the last complete round whatever its state,
#'   which is useful for diagnosing a diverging fit.
#'
#' @return A named list of (co)variance matrices, ending in \code{residuals},
#'   with an integer attribute \code{round}.
#'
#' @details
#' The result is not a drop-in for the \code{var.ini} argument of
#' \code{\link{remlf90}}, which covers only the terms of \code{random} plus the
#' residual; the genetic and spatial components have to go to their own slots:
#'
#' \preformatted{
#' vc <- reml_checkpoint("met3.log", model = fit)
#' remlf90(..., var.ini = vc[c("bl", "residuals")],
#'              genetic = list(..., var.ini = vc$genetic),
#'              spatial = list(..., var.ini = vc$spatial))
#' }
#'
#' To resume the \emph{same} model there is no need for any of this: use
#' \code{remlf90(..., progress_file = "met3.log", cont = TRUE)}.
#'
#' @seealso \code{\link{remlf90}}
#' @export
#' @examples
#' \dontrun{
#' ## From a second R session, while a fit is running:
#' vc <- reml_checkpoint("met3.log")
#' attr(vc, "round")
#' }
reml_checkpoint <- function(file, model = NULL, spd = TRUE) {

  if (!is.character(file) || length(file) != 1L || is.na(file))
    stop("'file' must be a single file name.", call. = FALSE)
  if (!file.exists(file))
    stop("'file' does not exist: ", file, call. = FALSE)

  x <- readLines(file, warn = FALSE)

  if (is.null(model)) {
    ## Without a model we have to infer how many groups to expect. Take the
    ## largest count any round shows, not the last one's: a round cut short by
    ## a crash shows too few, and believing it would quietly return a
    ## checkpoint with a component missing.
    anchors <- reml_round_index(x)
    if (!length(anchors))
      stop("'", file, "' holds no completed round.", call. = FALSE)

    grp_re <- "^[[:space:]]*(new [Gg]|G)[[:space:]]*$"
    bounds <- c(anchors, length(x) + 1L)
    n_groups <- max(vapply(
      seq_along(anchors),
      function(i) sum(grepl(grp_re,
                            x[seq.int(bounds[i], bounds[i + 1L] - 1L)])),
      1L))
    ck <- last_reml_checkpoint(x, n_groups = n_groups, spd = spd)
    if (!length(ck))
      stop("Cannot read a checkpoint from '", file, "': ",
           attr(ck, 'reason'), ".", call. = FALSE)
    return(ck)
  }

  if (!inherits(model, 'remlf90'))
    stop("'model' must be a fitted remlf90 object.", call. = FALSE)

  ans <- reml_checkpoint_var.ini(
    file, model$effects,
    ntraits = ncol(as.matrix(stats::model.response(model$mf))),
    lines = x, spd = spd)

  structure(ans$var, round = ans$round)
}


#' Validate the progress and resume arguments of remlf90()
#'
#' Runs before the binaries are checked so that the guards are exercisable
#' without the backend installed.
#'
#' @param progress_file character or NULL.
#' @param cont logical.
#' @param breedR.bin character. Location of the binaries, or the special values
#'   \code{"remote"} / \code{"submit"}.
#' @param debug logical. Resuming is refused under \code{debug}, which parses
#'   and returns nothing.
#'
#' @return the normalised \code{progress_file}, or \code{NULL}.
#' @keywords internal
check_progress_args <- function(progress_file, cont, breedR.bin,
                                debug = FALSE) {

  if (!is.logical(cont) || length(cont) != 1L || is.na(cont))
    stop("'cont' must be TRUE or FALSE.", call. = FALSE)

  ## Under debug the streaming branch still runs -- it is selected on
  ## progress_file, not on debug -- so a resume would rename the previous log
  ## aside, overwrite it, and then return NULL without parsing anything. There
  ## is nothing to recover from such a run, so refuse rather than consume the
  ## log for no result.
  if (isTRUE(cont) && isTRUE(debug))
    stop("'cont = TRUE' cannot be combined with 'debug = TRUE': a debug run ",
         "parses no results, so resuming would consume the previous log ",
         "without returning anything.", call. = FALSE)

  if (is.null(progress_file)) {
    if (isTRUE(cont))
      stop("'cont = TRUE' requires 'progress_file' pointing at the log of a ",
           "previous run.", call. = FALSE)
    return(NULL)
  }

  if (!is.character(progress_file) || length(progress_file) != 1L ||
      is.na(progress_file) || !nzchar(progress_file))
    stop("'progress_file' must be a single file name.", call. = FALSE)

  if (tolower(breedR.bin) %in% c("remote", "submit"))
    stop("'progress_file' and 'cont' require a local fit; they are not ",
         "available with breedR.bin = '", breedR.bin, "'.", call. = FALSE)

  ## remlf90() changes to a temporary directory before running the backend, so
  ## a relative path has to be resolved against the caller's working directory
  ## here or the log would silently be written into tempdir() and lost.
  dir <- dirname(progress_file)
  if (!dir.exists(dir))
    stop("The directory of 'progress_file' does not exist: ", dir,
         call. = FALSE)
  progress_file <- file.path(normalizePath(dir, winslash = "/",
                                           mustWork = TRUE),
                             basename(progress_file))

  if (isTRUE(cont) && !file.exists(progress_file))
    stop("'cont = TRUE' but 'progress_file' does not exist: ", progress_file,
         call. = FALSE)

  progress_file
}


#' Report a failed PROGSF90 run
#'
#' @param out character vector. Output of the backend.
#' @param status integer. Exit code.
#' @keywords internal
stop_progsf90_failure <- function(out, status) {
  stop("PROGSF90 binary failed with exit code ", status,
       ".\nOutput:\n", paste(utils::tail(out, 20), collapse = "\n"),
       call. = FALSE)
}
