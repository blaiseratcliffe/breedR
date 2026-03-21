## Tests for refactored code paths
## Verifies that behavioral changes and refactors produce correct results.

context("Refactored code paths")

## -- renderpf90.default stops on unsupported classes (L4) --

test_that("renderpf90 stops with a clear error on unsupported classes", {
  expect_error(renderpf90(42), "No renderpf90 method for class")
  expect_error(renderpf90("foo"), "No renderpf90 method for class")
  expect_error(renderpf90(list(a = 1)), "No renderpf90 method for class")
})

## -- renderpf90.blocks delegates to renderpf90.diagonal (M1) --

test_that("renderpf90.blocks produces identical output to renderpf90.diagonal", {
  ## Construct a minimal diagonal-like effect with a dgCMatrix incidence matrix
  n <- 10
  nlevels <- 3
  lvls <- sample(nlevels, n, replace = TRUE)
  mmx <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = lvls,
    x = rep(1, n),
    dims = c(n, nlevels)
  )

  ## Build objects that look like diagonal and blocks effects
  diag_effect <- list(incidence.matrix = mmx)
  class(diag_effect) <- "diagonal"

  block_effect <- list(incidence.matrix = mmx)
  class(block_effect) <- "blocks"

  result_diag <- renderpf90(diag_effect)
  result_block <- renderpf90(block_effect)

  expect_identical(result_diag, result_block)
})

## -- Vectorized build_grid indexing (Bolt) --

test_that("build_grid computes correct vector indices for a full grid", {
  coord <- expand.grid(row = 1:5, col = 1:4)
  grid <- build_grid(coord)

  ## Full 5x4 grid should produce sequential indices 1:20
  expect_equal(grid$idx, 1:20)
  expect_equivalent(grid$length, c(5, 4))
  expect_true(grid$regular)
})

test_that("build_grid computes correct indices for a sparse grid", {
  coord <- expand.grid(row = 1:5, col = 1:4)
  subset_idx <- c(1, 3, 7, 12, 20)
  coord2 <- coord[subset_idx, ]
  grid2 <- build_grid(coord2, autofill = FALSE)

  expect_equal(length(grid2$idx), 5)
  ## All indices should be within bounds of the grid
  expect_true(all(grid2$idx >= 1))
  expect_true(all(grid2$idx <= prod(grid2$length)))
  ## No duplicate indices for distinct coordinates
  expect_equal(length(unique(grid2$idx)), 5)
})

test_that("build_grid indices match column-major ordering", {
  ## For a 3x3 grid, column-major order is:
  ## (1,1)=1, (2,1)=2, (3,1)=3, (1,2)=4, (2,2)=5, ...
  coord <- expand.grid(row = 1:3, col = 1:3)
  grid <- build_grid(coord)

  ## Row 2, Col 3 should be index 2 + (3-1)*3 = 8
  row2col3 <- which(coord$row == 2 & coord$col == 3)
  expect_equal(grid$idx[row2col3], 8)
})

## -- Named constants are loaded and have expected values (L3) --

test_that("package constants have expected values", {
  expect_identical(MAX_COMPETITORS, 8L)
  expect_identical(MAX_REML_ITERATIONS, 5000L)
  expect_identical(SPLINE_BOUNDARY_KNOTS, 3L)
  expect_equal(REGULAR_GRID_QUANTILE_RANGE, c(.1, .6))
})
