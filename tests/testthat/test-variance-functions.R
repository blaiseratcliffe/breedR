## Tests for variance function helpers (R/varfunctions.R)

context("Variance function helpers")

## -- h2_formula() --

test_that("h2_formula produces correct simple heritability", {
  result <- h2_formula()
  expect_match(result, "^se_covar_function H2 ")
  expect_match(result, "G_2_2_1_1/\\(G_2_2_1_1\\+R_1_1\\)")
})

test_that("h2_formula includes other random effects in denominator", {
  result <- h2_formula(other_random = 3)
  expect_match(result, "G_2_2_1_1/\\(G_2_2_1_1\\+G_3_3_1_1\\+R_1_1\\)")
})

test_that("h2_formula handles multiple additional random effects", {
  result <- h2_formula(other_random = c(3, 4))
  expect_match(result, "G_3_3_1_1\\+G_4_4_1_1\\+R_1_1")
})

test_that("h2_formula handles maternal effect", {
  result <- h2_formula(maternal_effect = 3)
  # Denominator should include direct-maternal covariance and maternal variance
  expect_match(result, "G_2_3_1_1")
  expect_match(result, "G_3_3_1_1")
})

test_that("h2_formula handles multi-trait", {
  r1 <- h2_formula(trait = 1)
  r2 <- h2_formula(trait = 2)
  expect_match(r1, "G_2_2_1_1")
  expect_match(r2, "G_2_2_2_2")
  expect_match(r2, "R_2_2")
  expect_match(r2, "H2_t2")
})

test_that("h2_formula uses custom label", {
  result <- h2_formula(label = "my_h2")
  expect_match(result, "^se_covar_function my_h2 ")
})

test_that("h2_formula uses custom genetic effect number", {
  result <- h2_formula(genetic_effect = 5)
  expect_match(result, "G_5_5_1_1")
})

## -- rg_formula() --

test_that("rg_formula produces correct genetic correlation", {
  result <- rg_formula(1, 2)
  expect_match(result, "^se_covar_function rg_12 ")
  expect_match(result, "G_2_2_1_2/\\(G_2_2_1_1\\*G_2_2_2_2\\)\\*\\*0\\.5")
})

test_that("rg_formula handles different effect number", {
  result <- rg_formula(1, 3, effect = 4)
  expect_match(result, "G_4_4_1_3")
  expect_match(result, "G_4_4_1_1\\*G_4_4_3_3")
})

test_that("rg_formula uses custom label", {
  result <- rg_formula(1, 2, label = "gen_cor")
  expect_match(result, "^se_covar_function gen_cor ")
})

## -- vp_formula() --

test_that("vp_formula produces correct phenotypic variance", {
  result <- vp_formula()
  expect_match(result, "^se_covar_function Vp ")
  expect_match(result, "G_2_2_1_1\\+R_1_1")
})

test_that("vp_formula handles multiple random effects", {
  result <- vp_formula(random_effects = c(2, 3, 4))
  expect_match(result, "G_2_2_1_1\\+G_3_3_1_1\\+G_4_4_1_1\\+R_1_1")
})

test_that("vp_formula handles multi-trait", {
  result <- vp_formula(trait = 3)
  expect_match(result, "Vp_t3")
  expect_match(result, "G_2_2_3_3\\+R_3_3")
})

## -- var_ratio_formula() --

test_that("var_ratio_formula produces correct variance proportion", {
  result <- var_ratio_formula(effect = 3, all_random = c(2, 3))
  expect_match(result, "^se_covar_function prop_e3 ")
  expect_match(result, "G_3_3_1_1/\\(G_2_2_1_1\\+G_3_3_1_1\\+R_1_1\\)")
})

test_that("var_ratio_formula handles multi-trait", {
  result <- var_ratio_formula(effect = 3, trait = 2, all_random = c(2, 3))
  expect_match(result, "G_3_3_2_2")
  expect_match(result, "R_2_2")
})

## -- var_functions() --

test_that("var_functions generates h2 for single trait", {
  result <- var_functions(n_traits = 1)
  expect_length(result, 1)
  expect_match(result[1], "H2")
})

test_that("var_functions generates h2 + rg for two traits", {
  result <- var_functions(n_traits = 2)
  expect_length(result, 3)  # h2_t1, h2_t2, rg_12
  expect_match(result[1], "H2 ")
  expect_match(result[2], "H2_t2")
  expect_match(result[3], "rg_12")
})

test_that("var_functions generates correct count for three traits", {
  result <- var_functions(n_traits = 3, correlations = TRUE)
  # 3 heritabilities + 3 correlations (1-2, 1-3, 2-3) = 6
  expect_length(result, 6)
})

test_that("var_functions includes variance proportions for other random", {
  result <- var_functions(n_traits = 1, other_random = 3)
  expect_length(result, 2)  # h2 + prop_e3
  expect_match(result[2], "prop_e3")
})

test_that("var_functions skips correlations when requested", {
  result <- var_functions(n_traits = 2, correlations = FALSE)
  expect_length(result, 2)  # just 2 heritabilities
})

test_that("var_functions handles maternal + other random", {
  result <- var_functions(n_traits = 1, maternal_effect = 3, other_random = 4)
  # h2 (with maternal in denominator) + prop_e4
  expect_length(result, 2)
  expect_match(result[1], "G_2_3_1_1")  # maternal covariance in h2
})
