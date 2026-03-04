# test-de_silva_freq_EM.R
# Tests for the main de_silva_freq_EM() function.
# Gold standard: polysat::deSilvaFreq()

library(testthat)
library(polysat)
library(polyfreq)

# ---------------------------------------------------------------------------
# 1. Gold-standard comparison (from example_testing.R)
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM matches polysat::deSilvaFreq on simgen (self=0.1)", {
  skip_if_not_installed("polysat")
  result <- de_silva_freq_EM(simgen, self = 0.1)
  # Numeric frequencies must be identical to polysat within machine precision
  expect_equal(result, POLYSAT_RESULT, tolerance = 1e-8)
})

test_that("de_silva_freq_EM matches polysat::deSilvaFreq for self=0 (pure random mating)", {
  skip_if_not_installed("polysat")
  ps <- polysat::deSilvaFreq(simgen, self = 0)
  pf <- de_silva_freq_EM(simgen, self = 0)
  expect_equal(pf, ps, tolerance = 1e-8)
})

test_that("de_silva_freq_EM matches polysat::deSilvaFreq for self=1 (complete selfing)", {
  skip_if_not_installed("polysat")
  ps <- polysat::deSilvaFreq(simgen, self = 1)
  pf <- de_silva_freq_EM(simgen, self = 1)
  expect_equal(pf, ps, tolerance = 1e-8)
})

test_that("de_silva_freq_EM matches polysat::deSilvaFreq for self=0.5", {
  skip_if_not_installed("polysat")
  ps <- polysat::deSilvaFreq(simgen, self = 0.5)
  pf <- de_silva_freq_EM(simgen, self = 0.5)
  expect_equal(pf, ps, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 2. Output structure
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM returns a data frame", {
  result <- de_silva_freq_EM(simgen, self = 0.1)
  expect_s3_class(result, "data.frame")
})

test_that("de_silva_freq_EM output has correct row names (population names)", {
  result <- de_silva_freq_EM(simgen, self = 0.1)
  expect_equal(rownames(result), polysat::PopNames(simgen))
})

test_that("de_silva_freq_EM output has a Genomes column", {
  result <- de_silva_freq_EM(simgen, self = 0.1)
  expect_true("Genomes" %in% names(result))
})

test_that("de_silva_freq_EM output has null columns for each locus", {
  result <- de_silva_freq_EM(simgen, self = 0.1)
  null_cols <- paste0(polysat::Loci(simgen), ".null")
  expect_true(all(null_cols %in% names(result)))
})

test_that("allele frequencies sum to 1 per locus per population", {
  result   <- de_silva_freq_EM(simgen, self = 0.1)
  loci     <- polysat::Loci(simgen)
  for (L in loci) {
    cols <- grep(paste0("^", L, "\\."), names(result), value = TRUE)
    row_sums <- rowSums(result[, cols, drop = FALSE])
    expect_equal(row_sums, rep(1, nrow(result)), tolerance = 1e-8,
                 label = paste("Allele frequencies sum to 1 for locus", L))
  }
})

# ---------------------------------------------------------------------------
# 3. Input validation
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM errors when self is missing", {
  expect_error(de_silva_freq_EM(simgen), "Selfing rate required")
})

test_that("de_silva_freq_EM errors when ploidy is not uniform", {
  # Create a small genambig with mixed ploidies
  skip("Requires constructing a mixed-ploidy genambig — placeholder")
})

test_that("de_silva_freq_EM errors when ploidy is odd", {
  skip("Requires constructing an odd-ploidy genambig — placeholder")
})

test_that("de_silva_freq_EM errors when initFreq has no Genomes column", {
  bad_freq <- polysat::simpleFreq(simgen)
  bad_freq$Genomes <- NULL
  expect_error(de_silva_freq_EM(simgen, self = 0.1, initFreq = bad_freq),
               "Genomes")
})

# ---------------------------------------------------------------------------
# 4. Subset of samples / loci
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM accepts a samples subset", {
  sub_samples <- polysat::Samples(simgen, populations = "Pop1")
  ps <- polysat::deSilvaFreq(simgen, self = 0.1, samples = sub_samples)
  pf <- de_silva_freq_EM(simgen, self = 0.1, samples = sub_samples)
  expect_equal(pf, ps, tolerance = 1e-8)
})

test_that("de_silva_freq_EM accepts a single-locus subset", {
  ps <- polysat::deSilvaFreq(simgen, self = 0.1, loci = "loc1")
  pf <- de_silva_freq_EM(simgen, self = 0.1, loci = "loc1")
  expect_equal(pf, ps, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 5. initNull parameter
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM accepts a per-locus initNull vector", {
  loci     <- polysat::Loci(simgen)
  initNull <- setNames(c(0.1, 0.2, 0.05), loci)
  ps <- polysat::deSilvaFreq(simgen, self = 0.1, initNull = initNull)
  pf <- de_silva_freq_EM(simgen, self = 0.1, initNull = initNull)
  expect_equal(pf, ps, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 6. Convergence / iteration limits
# ---------------------------------------------------------------------------

test_that("de_silva_freq_EM respects maxiter", {
  # With maxiter=1 the result should differ from the converged solution
  ps_conv <- de_silva_freq_EM(simgen, self = 0.1)
  ps_1    <- de_silva_freq_EM(simgen, self = 0.1, maxiter = 1)
  expect_false(isTRUE(all.equal(ps_conv, ps_1, tolerance = 1e-3)))
})
