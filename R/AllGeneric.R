#' Covariance structure of a breedR component
#' 
#' This generic function returns the covariance or precision matrix of a breedR 
#' random effect or a group of effects.
#' 
#' For \code{effect_group}s, it returns the common structure of all the elements
#' in the group.
#' 
#' @param x A \code{breedr_effect}.
#'   
get_structure <- function(x) UseMethod('get_structure')

#' Parameters of a breedR component
#' 
#' @param x Some \code{breedR} component.
#'   
get_param <- function(x) UseMethod('get_param')

#' Type of a (group of) effect(s)
#' 
#' Generic function that returns whether an effect or a group of effects in a
#' breedR mixed model is \code{fixed} or \code{random}
#' @param x object to be \emph{translated} to progsf90
effect_type <- function(x) UseMethod('effect_type')

#' Render a progsf90 effect
#'
#' Translates breedR effects into progsf90 parameters and data.
#'
#' This is an internal function. Not exported.
#'
#' Each method must return a named list. The required fields depend on the
#' effect type:
#'
#' \strong{All effects} (fixed and random):
#' \describe{
#'   \item{pos}{integer vector. Column position(s) of this effect in the
#'     PROGSF90 data file. Length 1 for simple effects, length \code{n} for
#'     effects with \code{n} virtual sub-effects (e.g. competition has 8).}
#'   \item{levels}{integer vector, same length as \code{pos}. Number of levels
#'     for each virtual effect. For nested sub-effects, all but the last are 0.}
#'   \item{type}{character vector, same length as \code{pos}. Either
#'     \code{"cross"} (factor/indicator) or \code{"cov"} (covariate).}
#'   \item{nest}{integer vector or \code{NA}. Column positions of the nesting
#'     variables. \code{NA} when there is no nesting.}
#'   \item{data}{numeric matrix. The data columns to write to the PROGSF90
#'     data file for this effect.}
#' }
#'
#' \strong{Random effects} (diagonal, generic, ar, splines, blocks, genetic,
#' competition, pec) additionally require:
#' \describe{
#'   \item{model}{character string. The PROGSF90 model name, e.g.
#'     \code{"diagonal"}, \code{"add_animal"}, \code{"user_file"},
#'     \code{"user_file_i"}.}
#'   \item{file_name}{character string. Base name for the auxiliary structure
#'     file, or \code{""} when no file is needed (e.g. diagonal).}
#'   \item{file}{data.frame, matrix, or \code{""}. The content of the auxiliary
#'     file (pedigree data.frame, structure matrix in triplet format, etc.).}
#' }
#'
#' \strong{Effect groups} (\code{effect_group}) additionally include:
#' \describe{
#'   \item{var}{numeric matrix. The initial (co)variance matrix for the group.}
#' }
#'
#' @param x object of class breedr_modelframe, effect_group or breedr_effect.
#' @return A named list as described above.
#' @family renderpf90
renderpf90 <- function(x) UseMethod('renderpf90')

#' Get the Pedigree from an object
#' 
#' Returns an object from the formal class \code{pedigree}.
#' @param x object to extract pedigree from
#' @param ... Arguments to be passed to methods.
#' @references \code{\link[pedigreemm]{pedigree-class}} from package
#'   \code{pedigreemm}
#' @export
get_pedigree <- function(x, ...) UseMethod('get_pedigree')

#' Extract the number of traits
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
get_ntraits <- function(x, ...) UseMethod('get_ntraits')


#' 'move' an arrangement in a given direction
#' @param x matrix or list of matrices
#' @param dir a \emph{direction} in ('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')
neighbours.at <- function(x, dir) UseMethod('neighbours.at')


#' Number of generations
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
ngenerations <- function(x, ...) UseMethod('ngenerations')

#' Number of individuals
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
nindividuals <- function(x, ...) UseMethod('nindividuals')

#' Simulate a spatial structure
#' 
#' Takes a metagene simulated dataset and distributes its individuals randomly 
#' into a more or less square spatial region. Furthermore, it takes part of the 
#' phenotypic noise and puts some spatial structure into it.
#' 
#' Founders are not put into place, as they don't have phenotypic values. 
#' Therefore, they are returned without coordinates nor spatial values.
#' 
#' The variance of the spatial field is given as a proportion of the variance of
#' the random noise that was added to the Breeding Values to produce the 
#' phenotypes. The phenotypes are afterwards recalculated to preserve the 
#' heritability.
#' 
#' The spatial unit is the distance between consecutive trees.
#' @param meta A metagene object
#' @param variance A number between 0 and 1. The variance of the spatial field 
#'   as a proportion of the non-inheritable phenotypic variance. See Details.
#' @param range A number between 0 and 1. The range of the spatial field, as a 
#'   proportion of the region size.
#' @param ... Arguments to be passed to methods.
#' @return Another metagene object with spatial structure given through an 
#'   additional column \code{sp_X} with the spatially structured component of 
#'   the phenotype, and a 'spatial' list element with coordinates and simulation
#'   information
#' @export
sim.spatial <- function(meta, variance, range, ...) UseMethod('sim.spatial')
