
#### Context: breedR.setOption() ####
context("breedR options")

test_that("breedR.setOption() preserves a numeric vector value", {
  old <- breedR.getOption("ar.eval")
  on.exit(breedR.setOption("ar.eval", old), add = TRUE)
  new_grid <- c(-.9, -.3, .3, .9)
  breedR.setOption("ar.eval", new_grid)
  expect_equal(breedR.getOption("ar.eval"), new_grid)
})
