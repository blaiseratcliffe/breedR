#' Build an genetic model
#' 
#' Check conformity of arguments and return a \code{genetic} object.
#' 
#' This is a virtual class. No objects are expected to be created directly.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
genetic <- function(pedigree, incidence, covariance, precision) {
  
  ## checks
  stopifnot(!missing(pedigree))
  if (nrow(as.data.frame(pedigree)) != ncol(incidence))
    stop('The incidence matrix should have as many columns as individuals in the pedigree.')
  
  ## Build the random effect, and further specify the genetic class
  random.call <- mc <- match.call()
  arg.list <- as.list(mc)[-1]
  
  random.call[[1]] <- as.symbol('random')
  ans <- eval(random.call[-match('pedigree', as.list(random.call))],
              parent.frame())

  ans$pedigree <- pedigree
  class(ans) <- c('genetic', class(ans))
  
  return(ans)
}


#' Build an additive_genetic model
#' 
#' Check conformity of arguments and return a \code{additive_genetic} object.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' ped <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
#'                             dam  = c(NA,NA,2,NA,3,2),
#'                             label= 1:6)
#' inc <- cbind(0, 0, diag(4))
#' breedR:::additive_genetic(ped, inc)
additive_genetic <- function(pedigree, incidence) {
  
  
  ## Build the genetic effect, and further specify the additive_genetic class
  relationship.matrix <- pedigreemm::getA(pedigree)
  
  ## NOTE:
  ## We could as well use pedigreemm::getAInv(pedigree)
  ## to construct the precision matrix.
  
  ans <- genetic(pedigree, incidence, covariance = relationship.matrix)

  class(ans) <- c('additive_genetic', class(ans))
  
  return(ans)
}


## Internal codes of the individuals with original codes id
##
## Returns one integer per element of id, in order: the position of that
## individual in pedigree@label, which is its row in the pedigree and its
## column in a genetic incidence matrix. If build_pedigree() recoded the
## pedigree, id is in the original codes and is translated through
## attr(pedigree, 'map'). Stops, naming them, if any id is not the code of an
## individual in the pedigree: absent, 0, negative, fractional, NA or not
## numeric.
##
## Only a whole id within the map's range is used as a subscript. `[` would
## drop a 0, take a negative id as an exclusion and truncate a fractional
## one, so map[id] could return fewer codes than ids, or the wrong ones,
## without an error (#62). The map is not inverted with match(), which would
## hash a vector as long as the largest original code.
pedigree_index <- function(pedigree, id) {
  idx <- rep(NA_integer_, length(id))
  if (is.numeric(id)) {
    if (is.null(map <- attr(pedigree, 'map'))) {
      idx <- match(id, as.integer(pedigree@label))
    } else {
      ok <- which(id >= 1 & id <= length(map) & id == trunc(id))
      idx[ok] <- map[id[ok]]
    }
  }
  if (anyNA(idx))
    stop(paste('The following individuals in id are',
               'not represented in the pedigree:\n',
               toString(id[is.na(idx)])), call. = FALSE)
  idx
}


#' Build an additive-genetic animal model
#' 
#' Given a pedigree, and an index vector of observations, build and 
#' \code{additive_genetic_animal} model.
#' 
#' \code{idx} must hold the index of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this
#' function will convert the codes internally.
#' 
#' @param idx integer vector of observed individuals (in the original
#'   codification)
#' @inheritParams additive_genetic
#' @importFrom methods as
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id = 1:4,
#'                   sire = c(11, 11, 2, 3),
#'                   dam  = c(12, NA, 1, 12))
#' ped <- build_pedigree(1:3, data = dat)
#' breedR:::additive_genetic_animal(ped, dat$id)
additive_genetic_animal <- function(pedigree, idx) {
  
  ## Checks
  stopifnot(is.numeric(idx))
  # Not necessarily: might be multiple observations
  # stopifnot(length(idx) < nrow(as.data.frame(pedigree)))
  
  ## Incidence matrix
  ## It is possible that the pedigree has been recoded/reordered
  ## In that case, we need to recode the data file id codes as well.
  ## pedigree_index() returns one valid index into pedigree@label per id,
  ## or stops naming the ids that are not in the pedigree.
  idx <- pedigree_index(pedigree, idx)
  
  ## The pedigree might potentially have further individuals
  ## to evaluate (either founders, or descendants).
  inc.mat <- as(
    Matrix::sparseMatrix(i = seq_along(idx),
                         j = idx,
                         x = 1,
                         dims = c(length(idx),
                                  nrow(as.data.frame(pedigree)))),
    'indMatrix')
  
  ans <- additive_genetic(pedigree, inc.mat)
  class(ans) <- c('additive_genetic_animal', class(ans))
  return(ans)
}



#' Build an additive-genetic competition model
#' 
#' Return incidence and structure for a \code{additive_genetic_competition}
#' model, given the pedigree, the spatial coordinates and codes of the
#' observations and the competition decay parameter.
#' 
#' \code{id} must hold the codes of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this 
#' function will handle the codes internally.
#' 
#' @param id integer vector of numeric codes for observed individuals
#' @inheritParams additive_genetic
#' @inheritParams competition
#' @inheritParams build_grid
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id   = 1:5,
#'                   sire = c(11, 11, 2, 3, 2),
#'                   dam  = c(12, NA, 1, 12, 1),
#'                   x    = c(rep(1:2, times = 2), 3),
#'                   y    = c(rep(1:2, each = 2), 3))
#' ped <- build_pedigree(1:3, data = dat)
#' breedR:::additive_genetic_competition(ped, coord = dat[, c('x', 'y')], dat$id, 2)
additive_genetic_competition <- function(pedigree,
                                         coordinates,
                                         id,
                                         decay,
                                         autofill = TRUE) {
  
  ## Checks
  stopifnot(is.numeric(id))
  stopifnot(length(id) == nrow(coordinates))
  # Not necessarily: might be multiple observations per genotype
  # stopifnot( (n <- length(id)) < (p <- nrow(as.data.frame(pedigree))))
  
  n <- length(id)    # n observations
  p <- nrow(as.data.frame(pedigree))   # n genotypes
  
  ## additive_genetic_competition inherits from additive_genetic
  ## and from competition simultaneously.
  ## As S3 classes do not support multiple inheritance, simulate it
  ## manually by constructing separate competition and additive_genetic
  ## objects with dummy covariance and incidence matrices respectively
  ## and then composing manually the random effect
  cov.dummy <- Matrix::Diagonal(n)
  comp.aux <- competition(coordinates = coordinates, 
                          covariance  = cov.dummy,
                          decay       = decay,
                          autofill    = autofill)
  
  inc.dummy <- Matrix::Diagonal(p)
  ag.aux <- additive_genetic(pedigree, inc.dummy)
  
  ## the internal codes of the observed individuals are the indices
  ## of the corresponding levels of the random effect
  id.internal <- pedigree_index(pedigree, id)
  ## Competition incidence W Z: W (n x n) holds the weight of each neighbouring
  ## record, Z (n x p) maps each record to its individual, with zero columns for
  ## unobserved individuals (e.g. founders). The weights of several neighbouring
  ## records of the same individual (e.g. ramets of a clone) add up (#107).
  Z <- Matrix::sparseMatrix(i = seq_len(n), j = id.internal, x = 1, dims = c(n, p))
  inc.mat <- comp.aux$incidence.matrix %*% Z
  
  random.args <- structure(list(inc.mat, ag.aux$structure.matrix),
                           names = c('incidence', ag.aux$structure.type))
  
  ans <- do.call('random', random.args)
  ans$pedigree <- pedigree
  ans$coordinates <- coordinates
  
  ## Include inherited classes from the additive_genetic object 
  ## in the first place
  class(ans) <- c('additive_genetic_competition',
                  setdiff(class(ag.aux), class(ans)),
                  setdiff(class(comp.aux), class(ans)),
                  class(ans))
  
  return(ans)
}

