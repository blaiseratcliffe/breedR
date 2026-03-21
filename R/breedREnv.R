## Nothing to export

## Package-level constants
## Maximum number of competitors (8 cardinal/intercardinal directions on a grid)
MAX_COMPETITORS <- 8L
## Maximum REML iterations (hardcoded in the PROGSF90 Fortran binaries)
MAX_REML_ITERATIONS <- 5000L
## Number of boundary knots added at each end of B-spline basis
SPLINE_BOUNDARY_KNOTS <- 3L
## Quantile range used to detect regular grid spacing in fill_holes()
REGULAR_GRID_QUANTILE_RANGE <- c(.1, .6)

## Define the environment for breedR used to store options and the
## model-list. Reuse the environment if it is there already.

## Thanks to the INLA team, from where I took the whole option management system.

if (exists(".breedREnv") && is.environment(.breedREnv)) {
    ## then reuse it
} else {
    .breedREnv = new.env()
}

`breedR.get.breedREnv` = function(...)
{
    if (exists(".breedREnv") && is.environment(.breedREnv))
        return (.breedREnv)
    stop("Environment '.breedREnv' does not exists and is required for breedR to work. Restart 'R'.")
}

