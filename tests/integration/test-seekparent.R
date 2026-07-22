### Integration test: SNP-based parentage verification (requires binaries) ###

context("Parentage verification (seekparentf90)")

## Build a small Mendelian-consistent pedigree so conflicts are meaningful:
## 20 founders with random genotypes, 30 offspring inheriting one allele from
## each parent. One offspring is then given a deliberately wrong sire.
set.seed(5)
sk_nf <- 20L; sk_no <- 30L; sk_nsnp <- 400L
sk_fg <- matrix(sample(0:2, sk_nf * sk_nsnp, replace = TRUE,
                       prob = c(0.25, 0.5, 0.25)), nrow = sk_nf)
sk_sire <- sample(1:10, sk_no, replace = TRUE)
sk_dam  <- sample(11:20, sk_no, replace = TRUE)
sk_transmit <- function(g)
  ifelse(g == 0, 0, ifelse(g == 2, 1, sample(0:1, length(g), replace = TRUE)))
sk_og <- t(vapply(seq_len(sk_no),
                  function(i) sk_transmit(sk_fg[sk_sire[i], ]) +
                              sk_transmit(sk_fg[sk_dam[i], ]),
                  numeric(sk_nsnp)))
sk_G <- rbind(sk_fg, sk_og)
sk_ids <- c(seq_len(sk_nf), sk_nf + seq_len(sk_no))

sk_ped <- data.frame(animal = sk_ids,
                     sire = c(rep(0L, sk_nf), sk_sire),
                     dam  = c(rep(0L, sk_nf), sk_dam))
sk_wrong <- sk_nf + 1L
sk_ped$sire[sk_ped$animal == sk_wrong] <- 10L   # deliberately wrong sire

sk_snp <- file.path(tempdir(), "test_seek_snp.txt")
write_snp_file(sk_G, ids = sk_ids, file = sk_snp)
sk_pedf <- file.path(tempdir(), "test_seek_ped.txt")
write.table(sk_ped, sk_pedf, row.names = FALSE, col.names = FALSE, quote = FALSE)

sk_dir <- file.path(tempdir(), "test_seek_run")
unlink(sk_dir, recursive = TRUE)
sk_res <- seekparentf90(snp_file = sk_snp, ped_file = sk_pedf, dir = sk_dir)

test_that("seekparentf90 runs and returns conflict results", {
  expect_type(sk_res, "list")
  expect_false(is.null(sk_res$check))              # Check_<pedfile>
  expect_false(is.null(sk_res$conflicts))          # Parent_Progeny_Conflicts.txt
  expect_false(is.null(sk_res$conflicts_summary))  # ..._Summary.txt
})

test_that("seekparentf90 detects the injected wrong sire", {
  # The program reports the count of sire-progeny conflicts in its output.
  hit <- grep("Sire-progeny with conflicts", sk_res$output, value = TRUE)
  expect_length(hit, 1L)
  n_conf <- as.integer(regmatches(hit, regexpr("[0-9]+", hit)))
  expect_equal(n_conf, 1L)     # exactly the one animal we corrupted
  dam_hit <- grep("Dam-progeny with conflicts", sk_res$output, value = TRUE)
  expect_equal(as.integer(regmatches(dam_hit, regexpr("[0-9]+", dam_hit))), 0L)
})

test_that("seekparentf90 errors on a missing candidate-list file", {
  expect_error(
    seekparentf90(snp_file = sk_snp, ped_file = sk_pedf,
                  seek_sire = "no_such_candidates.txt",
                  dir = file.path(tempdir(), "test_seek_bad")),
    "seek_sire file not found")
})
