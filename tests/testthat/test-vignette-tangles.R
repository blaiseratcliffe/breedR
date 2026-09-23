## Regression tests for #30: inst/doc/Genomic-selection.R was shipped as an
## inert tangle (only the setup chunk live, all 7 illustrative chunks
## commented out) because it was produced by knitr's purl hook during a real
## rmarkdown::render(), where the vignette's own eval = FALSE default leaked
## into every later chunk.

test_that("the shipped Genomic-selection.R tangle has live code", {
  tangle_path <- system.file("doc", "Genomic-selection.R", package = "breedR")
  skip_if(identical(tangle_path, ""),
          "inst/doc/Genomic-selection.R not found via system.file()")

  lines <- readLines(tangle_path, warn = FALSE)
  live <- lines[!grepl("^\\s*(#.*)?$", lines)]

  expect_gt(length(live), 10)
  expect_true(any(grepl("^\\s*res\\s*<-\\s*remlf90\\(", live)))
  expect_false(any(grepl("opts_chunk\\$set\\([^)]*eval\\s*=\\s*FALSE", live)))
})

test_that("breedR.tangle_vignette() purls a vignette without leaking eval = FALSE", {
  skip_if_not(exists("breedR.tangle_vignette", mode = "function"),
              "breedR.tangle_vignette() not defined (only available via load_all(), not the installed package)")

  src <- tempfile(fileext = ".Rmd")
  out <- tempfile(fileext = ".R")
  on.exit(unlink(c(src, out)))

  writeLines(c(
    "---",
    "title: scratch",
    "---",
    "",
    "```{r setup, include = FALSE}",
    "knitr::opts_chunk$set(eval = FALSE)",
    "```",
    "",
    "```{r arith}",
    "1 + 1",
    "```"
  ), src)

  breedR.tangle_vignette(src, out)

  tangled <- readLines(out, warn = FALSE)
  expect_true(any(grepl("^1 \\+ 1$", tangled)))
  expect_false(any(grepl("^#\\s*1 \\+ 1$", tangled)))
})
