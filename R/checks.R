## Functions for checking model components
## Internal - not exported


#' Check initial variances specification
#' 
#' If the user specified initial values, verify that all random effects were
#' included. Otherwise, set default values. In any case, validate all initial
#' values.
#' 
#' @return A list with initial covariance matrices for all random effects in the
#'   model. A logical attribute `var.ini.default` is TRUE if values were set by
#'   default.
#' 
#' @param x list. user specification of var.ini (or NULL)
#' @param random formula. user specification of random effects.
#' @param response numeric vector or matrix.
#' @param traits internal named list of logical trait-presence masks.
#' @param from_checkpoint logical. Whether \code{remlf90()} is resuming from a
#'   progress file, which replaces every initial variance: defaults are then
#'   placeholders.
#' @return  matrix of observation values.
check_var.ini <- function (x, random, response, traits = NULL,
                           from_checkpoint = FALSE) {
  
  
  ## terms in the random component + 'residual'
  random.terms <- switch( is.null(random) + 1,
                          c(attr(stats::terms(random), 'term.labels'), 'residuals'),
                          'residuals')
  
  if (!is.null(x)) {
    
    ## set flag: initial variances not specified by user
    attr(x, 'var.ini.default') <- FALSE
    
    ## normalize names
    names(x) <- match.arg(names(x),
                          random.terms,
                          several.ok = TRUE)
    
    ## check that all the required variances are given
    if (!setequal(names(x), random.terms)) {
      stop(paste('Some initial variances missing in var.ini.\n',
                 'Please specify either all or none.'))
    }
  } else {
    
    ## set default values and flag
    div_fun <- breedR.getOption('default.initial.variance')
    default_ini <- if (isTRUE(from_checkpoint))
      resume_placeholder_variance(response)
      else eval(div_fun)(response, dim = 1, cor.effect = 0.1, digits = 2)
    x <- lapply(random.terms, function(x) default_ini)
    names(x) <- random.terms
    attr(x, 'var.ini.default') <- TRUE
  }
  
  ## validate values
  for (i in seq_along(x)) {
    active <- if (names(x)[i] != 'residuals') traits[[names(x)[i]]]
    if (is.null(active)) {
      validate_variance(x[[i]], what = names(x)[i], where = "var.ini specification")
    } else {
      validate_variance(x[[i]], dimension = rep(ncol(as.matrix(response)), 2),
                        what = names(x)[i], where = "var.ini specification",
                        active = active)
      x[[i]] <- mask_variance(x[[i]], active)
    }
  }
  
  ## return component with names normalised and 
  ## possibly default values added
  return(x)
}


## Checks and completes the specification of a genetic model
check_genetic <- function(model = c('add_animal', 'competition'),
                          pedigree,
                          id,
                          coordinates,
                          competition_decay = 1,
                          pec = FALSE,
                          autofill = TRUE,
                          var.ini,
                          data,
                          response,
                          trait.active = NULL,
                          pec.trait.active = NULL,
                          ...,
                          from_checkpoint = FALSE) {
  
  ## do not include data in the call
  ## data is an auxiliar for checking and substituting id
  ## but it is not part of the genetic component specification
  mc <- match.call()
  mc <- mc[!names(mc) %in% c('data', 'response', 'trait.active',
                             'pec.trait.active', 'from_checkpoint')]
  
  ## Mandatory arguments
  for (arg in c('model', 'pedigree', 'id')) {
    if (eval(call('missing', as.name(arg))) || 
        eval(call('is.null', as.name(arg))))
      stop(paste('Argument', arg, 'required in the genetic component.'))
  }
  
  ## Match argument model
  mc$model <- match.arg(model)
  
  ## Check type of argument pedigree 
  ## recode if necessary and always return a 'pedigree'
  if (!inherits(pedigree, 'pedigree')){
    ped.df <- try(as.data.frame(pedigree))
    if (inherits(ped.df, 'try-error')) 
      stop(paste('The argument pedigree in the genetic component',
                 'must be coercible to data.frame'))
    if (!ncol(ped.df) == 3)
      stop(paste('The argument pedigree in the genetic component',
                 'must have exactly 3 columns\ncorresponding to the',
                 'individual, its father and its mother respectively.'))
    pedigree <- build_pedigree(1:3, data = ped.df)
  }
  if (!all(check_pedigree(pedigree))) {
    pedigree <- build_pedigree(1:3, data = as.data.frame(pedigree))
  }
  mc$pedigree <- eval(pedigree, parent.frame())
  
  
  ## id must be either a variable name in data
  ## or a vector of codes in the pedigree
  if (length(id) == 1) {
    if (is.character(id) && id %in% names(data))
      mc$id <- as.integer(data[, id])
    else
      stop(paste('The argument id in the genetic component',
                 'must be either a vector of codes or a variable',
                 'name in the argument data.'))
  } else {
    mc$id <- id
  }

  ## The codes in id must correspond to valid codes in the pedigree
  ## possibly recoded
  if (!is.null(attr(mc$pedigree, "map")))
    recoded_id <- attr(mc$pedigree, "map")[mc$id]
  else recoded_id <- mc$id
  if (!all(idx <- recoded_id %in% mc$pedigree@label))
    stop(paste('The following individuals in id are',
               'not represented in the pedigree:\n',
               toString(mc$id[which(!idx)])))

  ## flag indicating whether the var.ini was taken by default
  ## or specified by the user
  attr(mc, 'var.ini.default') <- FALSE

  ## default initial variance function
  div_fun <- breedR.getOption('default.initial.variance')
  
  ## dimension of the genetic effect
  dim <- switch(mc$model, add_animal = 1, competition = 2)
  
  ## Set default var.ini if missing
  if (missing(var.ini) || is.null(var.ini)) {

    ## default initial covariance matrix
    var.ini <- if (isTRUE(from_checkpoint))
      resume_placeholder_variance(response, dim)
      else eval(div_fun)(response, dim = dim, cor.effect = 0.1, digits = 2)
    
    ## set flag indicating a default initial value
    attr(mc, 'var.ini.default') <- TRUE
  }
  
  ## Validate initial variance (SPD, dimensions, etc.)
  active <- if (!is.null(trait.active))
    expand_trait_mask(trait.active, ncol(as.matrix(response)), dim)
  validate_variance(
    var.ini,
    dimension = rep(dim*ncol(as.matrix(response)), 2),
    where = 'genetic component.', active = active)
  var.ini <- mask_variance(var.ini, active)
  
  ## Checks specific to competition models
  if (mc$model == 'competition') {
    
    ## Mandatory arguments
    for (arg in c('coordinates')) {
      if (eval(call('missing', as.name(arg))))
        stop(paste('Argument', arg, 'required in the genetic component.'))
    }
    
    mc$coordinates <- normalise_coordinates(coordinates,
                                            where = 'genetic component')
    
    ## Check pec argument
    # Specification of Permanent Environmental Effect
    ## Must be a list or a logical
    ## in the latter case, make it a list
    if (!is.list(pec)) {
      if (length(pec) == 1) {
        if (is.logical(pec)) {
          pec <- list(present = pec)
        } else {
          if (is.numeric(pec) && pec > 0) {
            pec <- list(present = TRUE, var.ini = pec)
          } else {
            stop('pec must be either list, a logical value or a positive number')
          }
        }
      } else {
        stop('pec must be either a list, a logical value or a positive number')
      }
    }
    
    ## Must be named
    if (is.null(names(pec)) || !all(nchar(names(pec))>0))
      stop('pec must be a named list')
    
    ## Match names
    names(pec) <- match.arg(names(pec), c('present', 'var.ini'), several.ok = TRUE)
    
    ## If there is no specification of 'present' it means it is present
    if (!'present' %in% names(pec)) {
      pec$present <- TRUE
    }
    
    ## Default initial variance
    if (!'var.ini' %in% names(pec)) {
      if (!attr(mc, 'var.ini.default') && pec$present) {
        stop(paste0('var.ini must be specified for pec as well, ',
                    'in the competition specification.\n',
                    'e.g. pec = list(present = TRUE, var.ini = 1)'))
      }
      
      ## default initial covariance matrix
      pec$var.ini <- if (isTRUE(from_checkpoint))
        resume_placeholder_variance(response)
        else eval(div_fun)(response, dim = 1, cor.effect = 0.1, digits = 2)
    }
    
    ## Validate initial variance in pec
    validate_variance(pec$var.ini,
                      what = "pec$var.ini",
                      where = "genetic component",
                      dimension = if (is.null(pec.trait.active)) dim(as.matrix(pec$var.ini))
                                  else rep(ncol(as.matrix(response)), 2),
                      active = pec.trait.active)
    pec$var.ini <- mask_variance(pec$var.ini, pec.trait.active)
    
    ## At this point, names should match exactly those
    if (!all(idx <- names(pec) %in% c('present', 'var.ini'))) {
      bad.args <- names(pec)[!idx]
      stop(paste0('Unrecognized argument',
                  ifelse(length(bad.args) == 1, '', 's'), ' ',
                  paste(bad.args, collapse = ', '), ' in pec'))
    }
    if (!is.logical(pec$present) | length(pec$present) != 1)
      stop('one logical value expected in pec$present')
    mc$pec <- pec
    
    ## TODO: check here whether none or all var.ini were specified
    ## and return var.ini.default as an attribute
    
    ## If missing, assume the default value
    stopifnot(is.numeric(competition_decay))
    stopifnot(competition_decay > 0)
    mc$competition_decay <- competition_decay
  }
  
  mc$var.ini <- var.ini
  mc$autofill <- autofill
  
  return(structure(as.list(mc[-1]),
                   var.ini.default = attr(mc, 'var.ini.default')))
}



check_spatial <- function(model = c('splines', 'AR', 'blocks'),
                          coordinates,
                          id,
                          n.knots,
                          rho,
                          autofill = TRUE,
                          sparse   = TRUE,
                          var.ini,
                          data,
                          response,
                          trait.active = NULL,
                          from_checkpoint = FALSE) {

  ## do not include data in the call
  ## data is an auxiliar for checking and substituting id
  ## but it is not part of the genetic component specification
  mc <- match.call()
  mc <- mc[!names(mc) %in% c('data', 'response', 'trait.active',
                             'from_checkpoint')]
  
  for (arg in c('model', 'coordinates')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the spatial component.'))
  }
  
  mc$model <- match.arg(model)
  
  mc$coordinates <- normalise_coordinates(coordinates, 'spatial component')
  
  ## If blocks model, include the values of the relevant covariate
  if (model == "blocks") {
    
    ## id must be either a variable name in data
    ## or a vector of codes in the pedigree
    if (length(id) == 1) {
      if (is.character(id) && id %in% names(data))
        mc$id <- as.integer(data[, id])
      else
        stop(paste('The argument id in the block component',
                   'must be either a vector of codes or a variable',
                   'name in the argument data.'))
    }
    mc$id <- eval(mc$id)
    
    # Only factors make sense for blocks
    # If it is already a factor, it may have
    # unobserved levels. Otherwise, make it a factor.
    if( !is.factor(mc$id) )
      mc$id <- as.factor(mc$id)
    
  }
  
  ## checks for splines models
  if (mc$model == 'splines') {
    if (!missing(n.knots)) {
      ## If n.knots specified, check consistency
      if (!is.vector(n.knots) || length(n.knots) !=2 || !all(n.knots%%1==0))
        stop(paste('n.knots must be a vector of two integers'))
      mc$n.knots <- n.knots
    }
  }
  
  ## checks for AR models
  if (model == 'AR'){
    
    ## rho not specified: make it NA in both dimensions
    if (missing(rho) || is.null(rho)) rho <- matrix(c(NA, NA), 1, 2)
    
    if (any(is.na(rho))) {
      ## A bare length-two vector is one row (one candidate pair, possibly with
      ## an NA to search that dimension's default grid) -- not two rows of a
      ## single column, which is what as.data.frame() would otherwise make of
      ## it inside build.AR.rho.grid().
      if (is.null(dim(rho)) && is.atomic(rho)) {
        if (length(rho) != 2)
          stop('rho must contain exactly two components')
        rho <- matrix(rho, nrow = 1, ncol = 2)
      }
      rho.grid <- build.AR.rho.grid(rho)
    } else {
      ## Fully specified: keep a plain vector as-is (a single fit, not a grid --
      ## downstream code tests is.null(nrow(spatial$rho)) to tell the two apart).
      ## A matrix/data.frame is an actual grid; convert it to a data.frame so
      ## transform(spatial$rho, ...) in remlf90() dispatches on
      ## transform.data.frame rather than falling through to transform.default,
      ## which evaluates in the wrong frame.
      rho.grid <- if (is.null(dim(rho))) rho else as.data.frame(rho)
    }
    

    check_rho_values <- function(rho) {
     if (!all(vapply(rho, is.numeric, TRUE)))
      stop('Argument rho in the spatial component must be numeric')
      if (any(abs(rho)>=1))
        stop('rho must contain numbers strictly between -1 and 1')
      if (!is.vector(rho))
        stop('rho must be a vector')
      if (length(rho)!=2)
        stop('rho must contain exactly two components')
      
      return(invisible(TRUE))
    }
    
    if (is.null(nrow(rho.grid))) {
      ## i.e. if is not really a grid
      check_rho_values(rho.grid)
    } else {
      ## grid case
      apply(rho.grid, 1, check_rho_values)
    }
    
    mc$rho <- rho.grid
  }
  
  ## flag indicating whether the var.ini was taken by default
  ## or specified by the user
  attr(mc, 'var.ini.default') <- FALSE
  
  ## default initial variance function
  div_fun <- breedR.getOption('default.initial.variance')

  ## dimension of the spatial effect
  dim <- 1
  
  if (missing(var.ini) || is.null(var.ini)) {
    
    ## default initial covariance matrix
    var.ini <- if (isTRUE(from_checkpoint))
      resume_placeholder_variance(response, dim)
      else eval(div_fun)(response, dim, cor.effect = 0.1, digits = 2)
    
    ## set flag indicating a default initial value
    attr(mc, 'var.ini.default') <- TRUE
  } 

  ## Validate initial variance (SPD, dimensions, etc.)
  validate_variance(
    var.ini,
    dimension = rep(dim*ncol(as.matrix(response)), 2),
    where = 'spatial component.', active = trait.active
  )
  mc$var.ini <- mask_variance(var.ini, trait.active)
  
  ## evaluate remaining parameters
  mc$autofill <- autofill
  mc$sparse   <- sparse
  
  return(structure(as.list(mc[-1]),
                   var.ini.default = attr(mc, 'var.ini.default')))
}



check_generic <- function(x, response, traits = NULL, from_checkpoint = FALSE){
  
  mc <- match.call()
  
  if (missing(x)) return(NULL)
  
  ## check general specification
  if (!is.list(x) || is.null(names(x)))
    stop('The generic component must be a named list.', call. = FALSE)
  if (!all(nchar(names(x))>0))
    stop('All elements of the generic component must be named.', call. = FALSE)
  if (any(duplicated(names(x))))
    stop('Duplicated names in generic elements.', call. = FALSE)
  if (!all(idx <- sapply(x,is.list))) {
    nm <- names(x)[!idx]
    if (length(nm) > 1)
      msg <- paste("Elements", paste(nm, collapse = ", "),
                   "of the generic component must be lists.")
    else
      msg <- paste("Element", paste(nm, collapse = ", "),
                   "of the generic component must be a list.")
    stop(msg, call. = FALSE)
  }
  
  ## validate individual elements
  for (arg.idx in seq_along(x)){ 
    id <- paste("generic component", names(x)[arg.idx])
    result <- do.call(
      'validate_generic_element', 
      c(x[[arg.idx]],
        response = list(response),
        where = id,
        trait.active = list(traits[[generic_effect_names(x)[arg.idx]]]),
        from_checkpoint = from_checkpoint)
    )
    ## If valid, the original spec might have been completed
    ## with a default initial variance
    x[[arg.idx]] <- result
  }
  
  ## Check default var.ini values
  ## Either all specified or all by default
  var.ini.default <- vapply(x, attr, TRUE, 'var.ini.default')
  if (any(var.ini.default) && any(!var.ini.default)) {
    stop(paste('Some initial variances missing in the generic component.\n',
               'Please specify either all or none.'), call. = FALSE)
  }
  
  ## Merge individual attributes into the list object
  for (i in seq_along(x)) attr(x[[i]], 'var.ini.default') <- NULL
  attr(x, 'var.ini.default') <- any(var.ini.default)
  
  return(x)
}


validate_generic_element <- function(incidence, 
                                     covariance, 
                                     precision, 
                                     var.ini, 
                                     response,
                                     where,
                                     trait.active = NULL,
                                     from_checkpoint = FALSE) {
  
  mc <- match.call()
  mc <- mc[!names(mc) %in% c('response', 'where', 'trait.active',
                             'from_checkpoint')]
  
  for (arg in c('incidence')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the', where), call. = FALSE)
  }
  if (!xor(missing(covariance), missing(precision)))
    stop(paste('Exactly one argument between covariance',
               'and precision must be specified in the', where), call. = FALSE)
  
  if (missing(covariance)) {
    structure <- precision
    str.name <- 'precision'
  }
  else {
    structure <- covariance
    str.name <- 'covariance'
  }
  if(!is.matrix(incidence) && !inherits(incidence, 'Matrix'))
    stop(paste('Argument incidence must be of type matrix in the', where),
         call. = FALSE)
  if(!is.matrix(structure) && !inherits(structure, 'Matrix'))
    stop(paste(str.name, 'must be of type matrix in the', where), call. = FALSE)
  if(ncol(incidence) != nrow(structure))
    stop(paste('Non conformant incidence and', str.name, 'matrices in the', where),
         call. = FALSE)

  ## flag indicating whether the var.ini was taken by default
  ## or specified by the user
  attr(mc, 'var.ini.default') <- FALSE
  
  ## default initial variance function
  div_fun <- breedR.getOption('default.initial.variance')
  
  ## dimension of the generic effect
  dim <- 1
  
  if (missing(var.ini) || is.null(var.ini)) {
    ## If not specified, return function that gives the value
    ## in order to check later whether the value is default or specified
    var.ini <- if (isTRUE(from_checkpoint))
      resume_placeholder_variance(response, dim)
      else eval(div_fun)(response, dim, cor.effect = 0.1, digits = 2)
    
    ## set flag indicating a default initial value
    attr(mc, 'var.ini.default') <- TRUE
  }
  
  ## Validate initial variance 
  ## even if default: the user could have changed the default function
  validate_variance(
    var.ini,
    dimension = rep(dim*ncol(as.matrix(response)), 2),
    where = where, active = trait.active)
  
  mc$var.ini <- mask_variance(var.ini, trait.active)
  
  return(structure(as.list(mc[-1]),
                   var.ini.default = attr(mc, 'var.ini.default')))
}


#' Normalise coordinates specification
#' 
#' If checks succeed, returns a complete normalised specification.
#' 
#' @param x matrix-like object to be checked
#' @param where string. Model component where coordinates were specified. For 
#'   error messages only. E.g. \code{where = 'genetic component'}.
#'   
#' @return a two-column data.frame, with numeric values.
normalise_coordinates <- function (x, where = '') {

    ## Check coordinates and cast to data.frame
  coord <- try(as.data.frame(x))
  if (inherits(coord, 'try-error') || !is.data.frame(coord))
    stop(paste('Argument coordinates in the', where,
              'not coercible to a data.frame'))

  ## Recast to data.frame. If nrow(coord) == 1, vapply returns
  ## a named vector, not a data.frame. E.g.:
  ## is.data.frame(vapply(data.frame(x=1, y=2), as.numeric, rep(1, 1)))
  coord <- try(as.data.frame(vapply(coord,
                                    as.numeric,
                                    rep(1, nrow(coord)))))
  
  if (inherits(coord, 'try-error'))
    stop(paste('Argument coordinates in the', where, 'not numeric'))
  if (ncol(coord) != 2)
    stop(paste('Only two dimensions admitted for coordinates',
               'in the', where))
  return(coord)
}


#' Check properties for a covariance matrix
#'
#' @param x number or matrix.
#' @param dimension numeric vector with dimensions of the matrix
#' @param what string. What are we validating
#' @param where string. Model component where coordinates were specified. For 
#'   error messages only. E.g. \code{where = 'competition specification'}.
#' @param active optional logical vector selecting the covariance coordinates
#'   that must form a positive-definite principal submatrix. The full matrix
#'   must still be finite, symmetric, and of the specified dimensions.
#'
#' @return \code{TRUE} if all checks pass
validate_variance <- function (x, dimension = dim(as.matrix(x)),
                               what = 'var.ini', where = '', active = NULL) {

  stopifnot(
    is.numeric(x <- as.matrix(x)),
    is.numeric(dimension),
    length(dimension) == 2
  )
  
  if (nrow(x)!=ncol(x))
    stop(paste(what, "must be a square matrix in the", where), call. = FALSE)
  if (length(x) != prod(dimension))
    stop(paste(what, "must be a", paste(dimension, collapse = 'x'),
               "matrix in the", where), call. = FALSE)
  if (!is.null(active)) {
    active <- trait_mask(active, nrow(x))
    if (!all(is.finite(x)) || !isSymmetric(x, check.attributes = FALSE))
      stop(paste(what, 'must be finite and symmetric in the', where), call. = FALSE)
    retained <- x[active, active, drop = FALSE]
    if (!all(eigen(retained, symmetric = TRUE, only.values = TRUE)$values > 0))
      stop("The active covariance block for '", what, "' must be SPD",
           if (!is.null(names(active))) paste0(' (traits: ',
             paste(unique(names(active)[active]), collapse = ', '), ')'),
           '.', call. = FALSE)
    return(TRUE)
  }
  ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(x, check.attributes = FALSE) || !all( ev > 0 ))
    stop(paste(what, "must be a SPD matrix in the", where), call. = FALSE)
  
  return(TRUE)
}


## Default initial variance under remlf90(cont = TRUE). The checkpoint
## replaces every initial variance after build.effects(), so the
## data-based default would be discarded, and it is undefined (all NA)
## when no record observes every trait (#71). Any SPD matrix of the right
## dimension passes the validations it meets first, and the checkpoint's
## dimension check is keyed on it.
resume_placeholder_variance <- function(response, dim = 1) {
  diag(dim * ncol(as.matrix(response)))
}


## Trait presence is independent of starting values and observation missingness.
## Keep these helpers shared by all effects, backend translation and recovery.
trait_mask <- function(trait.active, ntraits) {
  if (is.null(trait.active)) return(rep(TRUE, ntraits))
  if (!is.logical(trait.active) || length(trait.active) != ntraits ||
      anyNA(trait.active) || !any(trait.active))
    stop('Invalid internal trait-presence mask.', call. = FALSE)
  trait.active
}

expand_trait_mask <- function(trait.active, ntraits, size = 1L) {
  rep(trait_mask(trait.active, ntraits), times = size)
}

effect_trait_mask <- function(effect, ntraits) {
  trait_mask(effect$trait.active, ntraits)
}

effect_covariance_mask <- function(effect, ntraits) {
  expand_trait_mask(effect$trait.active, ntraits, length(effect$effects))
}

has_trait_restrictions <- function(effects) {
  any(vapply(effects, function(x) !is.null(x$trait.active) &&
               any(!x$trait.active), TRUE))
}

solution_trait_masks <- function(effects, ntraits) {
  masks <- lapply(effects, function(g)
    rep(list(effect_trait_mask(g, ntraits)), max(1L, length(g$effects))))
  stats::setNames(unlist(masks, recursive = FALSE), get_efnames(effects))
}

mask_variance <- function(x, active) {
  if (is.null(active) || all(active)) return(x)
  x <- as.matrix(x)
  x[!active, ] <- 0
  x[, !active] <- 0
  x
}

generic_effect_names <- function(generic) {
  nm <- names(generic)
  special <- nm %in% c('genetic', 'spatial')
  nm[special] <- paste0('generic_', nm[special])
  nm
}

random_group_names <- function(mf, genetic, spatial, generic) {
  tt <- attr(attr(mf, 'terms'), 'term.types')
  c(names(tt)[tt == 'random'],
    if (!is.null(genetic)) 'genetic',
    if (!is.null(genetic) && identical(genetic$model, 'competition') &&
        isTRUE(genetic$pec$present)) 'pec',
    if (!is.null(spatial)) 'spatial', generic_effect_names(generic))
}

check_effect_traits <- function(traits, response) {
  if (is.null(traits) || (is.list(traits) && !length(traits))) return(NULL)
  nm <- names(traits)
  if (!is.list(traits) || is.null(nm) || anyNA(nm) ||
      any(!nzchar(nm)) || anyDuplicated(nm))
    stop("'traits' must be a list with unique, nonempty random-effect group names.",
         call. = FALSE)
  tr <- colnames(response)
  if (ncol(response) < 2L || is.null(tr) || anyNA(tr) ||
      any(!nzchar(tr)) || anyDuplicated(tr))
    stop("'traits' requires a multivariate response with unique, nonempty column names.",
         call. = FALSE)
  lapply(stats::setNames(nm, nm), function(n) {
    selected <- traits[[n]]
    if (!is.character(selected) || !length(selected) || anyNA(selected) ||
        any(!nzchar(selected)) || anyDuplicated(selected))
      stop("traits[['", n, "']] must contain unique, nonmissing response names ",
           'and select at least one trait.', call. = FALSE)
    unknown <- setdiff(selected, tr)
    if (length(unknown))
      stop("Unknown trait '", paste(unknown, collapse = "', '"),
           "' for random-effect group '", n, "'. Available traits: ",
           paste(tr, collapse = ', '), '.', call. = FALSE)
    stats::setNames(tr %in% selected, tr)
  })
}

check_trait_groups <- function(traits, mf, genetic, spatial, generic) {
  if (is.null(traits)) return(NULL)
  groups <- random_group_names(mf, genetic, spatial, generic)
  tt <- attr(attr(mf, 'terms'), 'term.types')
  all_names <- c(names(tt)[tt == 'fixed'], groups)
  for (nm in names(traits)) {
    if (!nm %in% groups || sum(all_names == nm) != 1L)
      stop("Unknown or ambiguous random-effect group '", nm,
           "' in 'traits'. Available groups: ", paste(unique(groups), collapse = ', '),
           '. For different fixed effects by trait, use renumf90().', call. = FALSE)
  }
  traits <- traits[!vapply(traits, all, TRUE)]
  if (!length(traits)) NULL else traits
}
