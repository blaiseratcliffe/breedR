### Remote execution plumbing ###

context("Remote computing")

## The remote path needs a configured Linux server, so nothing here fits a
## model. What can be checked without one is the property the fix relies on.

test_that("results are never retrieved into a shared directory", {

  ## retrieve_remote() used to default `dest` to tempdir(), and breedR.qget()
  ## took that default: every job of a session was untarred into the same
  ## place, so a second retrieval overwrote the first one's LOG and solutions
  ## and both objects then advertised the one path. Every caller now passes a
  ## directory of its own, and `dest` has no default so that a third caller
  ## cannot silently go back to sharing one.
  expect_true(identical(formals(retrieve_remote)$dest, quote(expr = )))
  expect_true(identical(formals(breedR.remote)$dest, quote(expr = )))

  ## and the one caller that has no per-fit directory to pass makes itself one
  expect_match(paste(deparse(body(breedR.qget)), collapse = ' '),
               'retrieve_remote\\(rdir, dest = breedR_workdir')
})
