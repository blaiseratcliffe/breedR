# breedR 0.13-0

* `remlf90()` accepts `traits`, a named list selecting the responses on which
  each random-effect group is fitted (#51). Unlisted groups retain all traits.
  Full covariance dimensions and existing fitted-model methods are preserved;
  absent coefficients and variances are `NA`, with zero contribution to fitted
  values. Checkpoints retain and validate the same selection. Fixed-effect
  selection and restrictions on the genetic group with `genomic` are not
  supported. Calls without a selection are unchanged.
