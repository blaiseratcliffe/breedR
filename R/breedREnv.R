## Nothing to export

## Package-level constants
## Maximum number of competitors (8 cardinal/intercardinal directions on a grid)
MAX_COMPETITORS <- 8L
## Default REML iteration cap of BLUPF90+ (OPTION maxrounds): a fit that
## cannot converge stops at round 10000, and 2.73/2.76 echo "default=10000".
## The 5000 of the manual and --help belongs to the legacy REMLF90/AIREMLF90.
MAX_REML_ITERATIONS <- 10000L
## Default REML convergence criterion of BLUPF90+ (OPTION conv_crit), as its
## output echoes it: "convergence criterion (default=1e-12)".
REML_CONV_CRIT <- 1e-12
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

