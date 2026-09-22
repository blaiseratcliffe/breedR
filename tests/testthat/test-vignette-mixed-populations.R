## Regression test for #36: method 3 of the Additive-Genetic vignette scaled the
## hybrid blocks of A by averages of the two population variances. Those factors
## are the covariance implied by the simulation only at lambda = 1: a hybrid and
## a pure-E relative covary by sigma2_E * A, not (3+lambda)/4 * sigma2_E * A.

context("Vignette: additive genetic models in mixed populations")

## Pull a function definition out of a tangled vignette, by brace matching
extract_definition <- function(path, name) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^", name, " <- function"), lines)
  if (length(start) != 1L) return(NULL)
  tail_lines <- lines[start:length(lines)]
  opens  <- vapply(gregexpr("[{]", tail_lines), function(m) sum(m > 0), 1L)
  closes <- vapply(gregexpr("[}]", tail_lines), function(m) sum(m > 0), 1L)
  end <- which(cumsum(opens - closes) == 0 & cumsum(opens) > 0)[1]
  tail_lines[seq_len(end)]
}

test_that("the vignette scale_A() gives the mixed-population covariance", {

  tangle <- system.file("doc", "Additive-Genetic-Models-in-Mixed-Populations.R",
                        package = "breedR")
  skip_if(identical(tangle, ""), "vignette tangle not found via system.file()")

  def <- extract_definition(tangle, "scale_A")
  expect_false(is.null(def))

  ## A small mixed population in the terms the vignette uses: 3 E founders, 2 J
  ## founders, and offspring covering every kind of pair -- pure, selfed, hybrid
  ## with pure, and hybrids sharing their E parent, their J parent, or both.
  n.founders <- c(E = 3, J = 2)
  founders <- data.frame(id = c("E1", "E2", "E3", "J1", "J2"),
                         pop.idx = c(1, 1, 1, 2, 2))
  crosses <- rbind(c(1, 2), c(1, 3), c(3, 3),             # EE, the last selfed
                   c(1, 4), c(1, 4), c(2, 4), c(1, 5),    # EJ
                   c(4, 5), c(4, 4))                      # JJ
  dat <- data.frame(
    id  = nrow(founders) + seq_len(nrow(crosses)),
    dad = crosses[, 1],
    mum = crosses[, 2],
    sp  = factor(apply(crosses, 1, function(x)
      paste(names(n.founders)[founders$pop.idx[sort(x)]], collapse = ""))))

  ped <- build_pedigree(1:3, data = dat)
  A <- pedigreemm::getA(ped)
  idx_pop <- function(x) {
    if (nchar(x) == 1) grep(x, founders$id)
    else match(dat$id[dat$sp == x], as.data.frame(ped)$self)
  }

  eval(parse(text = def))          # defines scale_A() against these objects

  ## Covariance implied by the simulation: bv = mean(parent BVs) + msp, with
  ## Var(msp) = (s2_dad + s2_mum)/4 and independent founders.
  P <- matrix(0, nrow(dat), nrow(founders))
  for (i in seq_len(nrow(dat))) {
    P[i, dat$dad[i]] <- P[i, dat$dad[i]] + .5
    P[i, dat$mum[i]] <- P[i, dat$mum[i]] + .5
  }
  Pf <- rbind(diag(nrow(founders)), P)
  Sigma <- function(s2) {
    s2f <- s2[founders$pop.idx]
    Pf %*% diag(s2f) %*% t(Pf) +
      diag(c(rep(0, nrow(founders)), (s2f[dat$dad] + s2f[dat$mum]) / 4))
  }

  for (lambda in c(0.3, 2/3, 1, 1.5)) {
    expect_equal(unname(as.matrix(3 * scale_A(lambda))),
                 unname(Sigma(c(3, 3 * lambda))),
                 tolerance = 1e-12,
                 info = paste("lambda =", lambda))
  }
})
