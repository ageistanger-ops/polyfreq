# test-helpers.R
# Unit tests for each internal helper function:
#   G_r, INDEXG_r, GENLIST_r, RANMUL_r, SELFMAT_r

library(testthat)
library(polyfreq)

# ---------------------------------------------------------------------------
# G_r  (must equal choose(n + q, q + 1))
# ---------------------------------------------------------------------------

test_that("G_r(0, n) == n for several values of n", {
  for (n in c(1L, 2L, 5L, 10L)) {
    expect_equal(G_r(0L, n), n, label = paste0("G_r(0, ", n, ")"))
  }
})

test_that("G_r(q, n) equals choose(n + q, q + 1)", {
  for (q in 0:4) {
    for (n in c(1L, 3L, 5L, 8L)) {
      expect_equal(G_r(q, n), as.integer(choose(n + q, q + 1)),
                   label = paste0("G_r(", q, ", ", n, ")"))
    }
  }
})

test_that("G_r matches polysat C++ G for tetraploid smatdiv input", {
  # smatdiv = (G(m-1, m+1))^2 where m = m2/2
  # For tetraploid (m2=4, m=2): G(1, 3) = choose(4, 2) = 6 -> smatdiv = 36
  expect_equal(G_r(1L, 3L), 6L)
  # For hexaploid (m2=6, m=3): G(2, 4) = choose(6, 3) = 20 -> smatdiv = 400
  expect_equal(G_r(2L, 4L), 20L)
  # For octoploid (m2=8, m=4): G(3, 5) = choose(8, 4) = 70 -> smatdiv = 4900
  expect_equal(G_r(3L, 5L), 70L)
})

# ---------------------------------------------------------------------------
# INDEXG_r
# ---------------------------------------------------------------------------

test_that("INDEXG_r returns 1 for the first (all-ones) genotype", {
  # In GENLIST, the first genotype is always (1,1,...,1)
  for (m2 in c(2L, 4L, 6L)) {
    expect_equal(INDEXG_r(rep(1L, m2), na1 = 3L, m2 = m2), 1L,
                 label = paste0("m2=", m2))
  }
})

test_that("INDEXG_r is consistent with GENLIST_r row ordering", {
  # Every row of GENLIST_r should map back to its own row index via INDEXG_r
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  for (g in seq_len(ng)) {
    expect_equal(INDEXG_r(ag[g, ], na1, m2), g,
                 label = paste0("row ", g))
  }
})

test_that("INDEXG_r is consistent with GENLIST_r for hexaploid", {
  na1 <- 4L; m2 <- 6L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  for (g in seq_len(ng)) {
    expect_equal(INDEXG_r(ag[g, ], na1, m2), g)
  }
})

# ---------------------------------------------------------------------------
# GENLIST_r
# ---------------------------------------------------------------------------

test_that("GENLIST_r returns a matrix with ng rows and m2 columns", {
  na1 <- 4L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  expect_equal(nrow(ag), ng)
  expect_equal(ncol(ag), m2)
})

test_that("GENLIST_r rows are non-decreasing (sorted genotypes)", {
  na1 <- 5L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  sorted <- apply(ag, 1, function(r) all(diff(r) >= 0))
  expect_true(all(sorted))
})

test_that("GENLIST_r allele indices are in 1:na1", {
  na1 <- 4L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  expect_true(all(ag >= 1L & ag <= na1))
})

test_that("GENLIST_r produces exactly ng = choose(na1+m2-1, m2) rows", {
  for (na1 in c(2L, 3L, 5L)) {
    for (m2 in c(2L, 4L, 6L)) {
      ng <- as.integer(choose(na1 + m2 - 1, m2))
      ag <- GENLIST_r(ng, na1, m2)
      expect_equal(nrow(ag), ng, label = paste0("na1=", na1, " m2=", m2))
    }
  }
})

test_that("GENLIST_r matches polysat GENLIST output (tetraploid, na1=3)", {
  skip_if_not_installed("polysat")
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  # polysat GENLIST via the exported C++ function (if available):
  if (existsMethod("GENLIST", "polysat")) {
    ps_ag <- polysat::GENLIST(ng, na1, m2)
    pf_ag <- GENLIST_r(ng, na1, m2)
    expect_equal(pf_ag, ps_ag)
  } else {
    skip("polysat GENLIST not directly accessible")
  }
})

# ---------------------------------------------------------------------------
# RANMUL_r
# ---------------------------------------------------------------------------

test_that("RANMUL_r returns a list with rmul and arep", {
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  expect_named(res, c("rmul", "arep"))
})

test_that("RANMUL_r rmul has length ng", {
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  expect_equal(length(res$rmul), ng)
})

test_that("RANMUL_r arep has dimensions ng x na1", {
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  expect_equal(dim(res$arep), c(ng, na1))
})

test_that("RANMUL_r arep row sums equal m2 (copy counts sum to ploidy)", {
  na1 <- 4L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  expect_equal(rowSums(res$arep), rep(m2, ng))
})

test_that("RANMUL_r rmul is m2! for all-distinct genotypes", {
  # A genotype with all distinct alleles has multinomial coeff = m2!
  # e.g. for m2=4, a quadriallelic genotype: rmul should be 24
  na1 <- 5L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  # Find quadriallelic genotypes (4 distinct alleles)
  distinct_count <- apply(ag, 1, function(r) length(unique(r)))
  quadri_rows    <- which(distinct_count == m2)
  expect_true(all(res$rmul[quadri_rows] == factorial(m2)))
})

test_that("RANMUL_r rmul is 1 for homozygous genotypes", {
  # All-same-allele genotype: multinomial coeff = m2! / m2! = 1
  na1 <- 4L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  res <- RANMUL_r(ng, na1, ag, m2)
  homo_rows <- which(apply(ag, 1, function(r) length(unique(r)) == 1))
  expect_true(all(res$rmul[homo_rows] == 1L))
})

test_that("RANMUL_r matches polysat RANMUL output", {
  skip_if_not_installed("polysat")
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  pf  <- RANMUL_r(ng, na1, ag, m2)
  if (exists("RANMUL", envir = asNamespace("polysat"))) {
    ps <- polysat:::RANMUL(ng, na1, ag, m2)
    expect_equal(pf$rmul, ps$rmul)
    expect_equal(pf$arep, ps$arep)
  } else {
    skip("polysat RANMUL not accessible")
  }
})

# ---------------------------------------------------------------------------
# SELFMAT_r
# ---------------------------------------------------------------------------

test_that("SELFMAT_r returns a square matrix of dimensions ng x ng", {
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  sm  <- SELFMAT_r(ng, na1, ag, m2)
  expect_equal(dim(sm), c(ng, ng))
})

test_that("SELFMAT_r row sums equal smatdiv = choose(m2, m)^2", {
  # Each row of the unnormalized selfing matrix must sum to the total number
  # of gamete pairs = C(m2, m)^2
  na1 <- 3L; m2 <- 4L; m <- m2 %/% 2L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  sm  <- SELFMAT_r(ng, na1, ag, m2)
  smatdiv <- choose(m2, m)^2
  expect_equal(rowSums(sm), rep(smatdiv, ng), tolerance = 0)
})

test_that("SELFMAT_r is non-negative integer", {
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  sm  <- SELFMAT_r(ng, na1, ag, m2)
  expect_true(all(sm >= 0L))
  expect_true(is.integer(sm) || is.numeric(sm))
})

test_that("SELFMAT_r matches polysat SELFMAT output (tetraploid, na1=3)", {
  skip_if_not_installed("polysat")
  na1 <- 3L; m2 <- 4L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  pf  <- SELFMAT_r(ng, na1, ag, m2)
  if (exists("SELFMAT", envir = asNamespace("polysat"))) {
    ps <- polysat:::SELFMAT(ng, na1, ag, m2)
    expect_equal(pf, ps)
  } else {
    skip("polysat SELFMAT not accessible")
  }
})

test_that("SELFMAT_r matches polysat SELFMAT output (hexaploid, na1=4)", {
  skip_if_not_installed("polysat")
  na1 <- 4L; m2 <- 6L
  ng  <- as.integer(choose(na1 + m2 - 1, m2))
  ag  <- GENLIST_r(ng, na1, m2)
  pf  <- SELFMAT_r(ng, na1, ag, m2)
  if (exists("SELFMAT", envir = asNamespace("polysat"))) {
    ps <- polysat:::SELFMAT(ng, na1, ag, m2)
    expect_equal(pf, ps)
  } else {
    skip("polysat SELFMAT not accessible")
  }
})
