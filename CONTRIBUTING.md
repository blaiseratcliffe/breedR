# Contributing to breedR

The most useful thing you can do is run breedR on your own data and compare the estimates with whatever you use now: BLUPF90 run directly, ASReml-R or sommer, for example. If anything diverges, [open an issue](https://github.com/blaiseratcliffe/breedR/issues) with the model call, both sets of estimates, and your R and OS versions.

## When the backend fails

If a BLUPF90 program fails, breedR's error message ends with the last 20 lines of its output. Paste those into the issue.

For the full log, rerun a local fit with `progress_file = "run.log"`. breedR then writes the backend's output to `run.log` as it is produced, and its error stream, if there is any, to `run.log.err`. Both files survive the error. `progress_file` needs a single model run, so if you are fitting an AR spatial model, fix `rho` to one pair first: the default `rho` grid search fits one model per value and refuses `progress_file`.

## Confidential data

A synthetic reproduction is just as useful as the real data. `breedR.sample.phenotype()` simulates phenotypes from a model and variances you specify (see `?breedR.sample.phenotype`), so you can often rebuild the structure of your data set without sharing any of it.

## Code changes

Pull requests are welcome. Run the unit tests with `Rscript -e "devtools::test()"` before opening one, and add a line at the top of `NEWS` for any change users will notice.
