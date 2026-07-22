### Tests for download-integrity hardening in binaries.R ###

context("Binary download integrity")

test_that("looks_like_executable accepts a valid PE header and rejects junk", {
  # A fake but structurally-valid Windows executable: "MZ" + padding > 1 KiB
  win <- tempfile(fileext = ".exe")
  writeBin(c(as.raw(c(0x4D, 0x5A)), as.raw(rep(0L, 2048))), win)
  expect_true(looks_like_executable(win, "windows"))
  expect_false(looks_like_executable(win, "linux"))   # wrong magic for ELF

  # A valid ELF header
  elf <- tempfile()
  writeBin(c(as.raw(c(0x7F, 0x45, 0x4C, 0x46)), as.raw(rep(0L, 2048))), elf)
  expect_true(looks_like_executable(elf, "linux"))

  # An HTML error page must be rejected
  html <- tempfile(fileext = ".exe")
  writeLines(rep("<html>404</html>", 100), html)
  expect_false(looks_like_executable(html, "windows"))

  # A file below the size floor is rejected even with the right magic
  tiny <- tempfile()
  writeBin(as.raw(c(0x4D, 0x5A)), tiny)
  expect_false(looks_like_executable(tiny, "windows"))
})

test_that("breedr_progsf90_repo validates the PROGSF90_URL scheme", {
  old <- Sys.getenv("PROGSF90_URL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("PROGSF90_URL")
          else Sys.setenv(PROGSF90_URL = old))

  Sys.setenv(PROGSF90_URL = "https://example.org/x")
  expect_equal(breedr_progsf90_repo(), "https://example.org/x")

  Sys.setenv(PROGSF90_URL = "/local/mirror")
  expect_equal(breedr_progsf90_repo(), "/local/mirror")

  Sys.setenv(PROGSF90_URL = "http://insecure.example/x")
  expect_warning(breedr_progsf90_repo(), "insecure")

  Sys.setenv(PROGSF90_URL = "ftp://evil.example/x")
  expect_error(breedr_progsf90_repo(), "https")

  Sys.unsetenv("PROGSF90_URL")
  expect_match(breedr_progsf90_repo(), "^https://")
})
