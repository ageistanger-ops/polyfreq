# helpers_desilva.R
#
# Pure-R implementations of the de Silva et al. (2005) subroutines.
# These replace the C/C++ compiled functions used in the polysat package
# (G, INDEXG, GENLIST, RANMUL, SELFMAT).
#
# Mathematical identity used throughout:
#   G(q, n)  =  choose(n + q, q + 1)
#
# Reference: De Silva HN et al. (2005) Heredity 95:327-334.
#            doi:10.1038/sj.hdy.6800728


# ---------------------------------------------------------------------------
# G_r
# ---------------------------------------------------------------------------

#' Generalized binomial coefficient (de Silva G function)
#'
#' Computes the value of the G function used throughout the de Silva et al.
#' (2005) algorithm.  Mathematically identical to \code{choose(n + q, q + 1)}.
#'
#' The original SAS/IML and C++ implementations compute this as a running
#' product:
#' \deqn{G(q, n) = \prod_{j=0}^{q} \frac{n + j}{j + 1}
#'               = \frac{n(n+1)\cdots(n+q)}{(q+1)!}
#'               = \binom{n+q}{q+1}}
#'
#' @param q Non-negative integer; upper summation index.
#' @param n Integer; base value.
#'
#' @return An integer scalar equal to \code{choose(n + q, q + 1)}.
#'
#' @examples
#' G_r(0, 5)  # 5
#' G_r(1, 5)  # 15
#' G_r(2, 5)  # 35
#'
#' @export
G_r <- function(q, n) {
  as.integer(choose(n + q, q + 1L))
}


# ---------------------------------------------------------------------------
# INDEXG_r
# ---------------------------------------------------------------------------

#' Map a genotype vector to its unique integer index
#'
#' Translates the INDEXG subroutine from de Silva et al. (2005).  Given a
#' sorted genotype vector \code{ag1} (allele indices in \code{1:na1}, length
#' \code{m2}), returns the position of that genotype in the lexicographic
#' enumeration produced by \code{\link{GENLIST_r}}.
#'
#' @param ag1 Integer vector of length \code{m2}; sorted allele indices
#'   (\eqn{1 \le ag1[i] \le na1}, non-decreasing).
#' @param na1 Integer; total number of alleles including the null
#'   (\code{na + 1}).
#' @param m2 Integer; ploidy level (must be even).
#'
#' @return An integer scalar giving the 1-based row index of genotype
#'   \code{ag1} in the matrix returned by \code{GENLIST_r}.
#'
#' @seealso \code{\link{GENLIST_r}}, \code{\link{G_r}}
#'
#' @examples
#' # Tetraploid (m2 = 4), 3 alleles + null (na1 = 4)
#' INDEXG_r(c(1L, 1L, 1L, 1L), na1 = 4L, m2 = 4L)  # 1
#' INDEXG_r(c(1L, 1L, 1L, 2L), na1 = 4L, m2 = 4L)  # 2
#'
#' @export
INDEXG_r <- function(ag1, na1, m2) {
  # C++ ag1[m2-1] -> R ag1[m2], C++ ag1[m2-2] -> R ag1[m2-1]
  x <- 1L + ag1[m2] - ag1[m2 - 1L]

  if (m2 >= 3L) {
    for (q in 1L:(m2 - 2L)) {
      # C++ ag1[m2-q-2] -> R ag1[m2-q-1], C++ ag1[m2-q-1] -> R ag1[m2-q]
      x <- x + G_r(q, na1 + 1L - ag1[m2 - q - 1L]) -
               G_r(q, na1 + 1L - ag1[m2 - q])
    }
  }

  # C++ ag1[0] -> R ag1[1]
  x <- x + G_r(m2 - 1L, na1) - G_r(m2 - 1L, na1 + 1L - ag1[1L])
  as.integer(x)
}


# ---------------------------------------------------------------------------
# GENLIST_r
# ---------------------------------------------------------------------------

#' Generate the full list of polyploid genotypes
#'
#' Translates the GENLIST subroutine from de Silva et al. (2005).  Enumerates
#' all \code{ng} possible genotypes for a polyploid of ploidy \code{m2} with
#' \code{na1} alleles (including null) as an integer matrix whose rows are
#' sorted allele-index vectors.
#'
#' Genotypes are combinations-with-replacement of \code{1:na1} taken \code{m2}
#' at a time, listed in lexicographic (non-decreasing) order.  The total
#' number of genotypes is:
#' \deqn{ng = \binom{na1 + m2 - 1}{m2} = G(m2 - 1,\, na1)}
#'
#' @param ng Integer; total number of genotypes (must equal
#'   \code{choose(na1 + m2 - 1, m2)}).
#' @param na1 Integer; number of alleles including the null allele.
#' @param m2 Integer; ploidy level.
#'
#' @return An integer matrix of dimensions \code{ng} x \code{m2}.  Each row
#'   is one genotype: a non-decreasing sequence of allele indices in
#'   \code{1:na1}.
#'
#' @seealso \code{\link{INDEXG_r}}, \code{\link{RANMUL_r}}
#'
#' @examples
#' # All tetraploid genotypes with 2 alleles + null (na1 = 3)
#' GENLIST_r(ng = 15L, na1 = 3L, m2 = 4L)
#'
#' @export
GENLIST_r <- function(ng, na1, m2) {
  ag  <- matrix(0L, nrow = ng, ncol = m2)
  ag1 <- rep(1L, m2)
  ag[1L, ] <- ag1

  g <- 1L
  a <- m2                       # 1-based (rightmost column)
  while (a > 0L) {
    if (ag1[a] == na1) {
      a <- a - 1L
    } else {
      ag1[a] <- ag1[a] + 1L
      if (a < m2) {
        ag1[(a + 1L):m2] <- ag1[a]
      }
      g <- g + 1L
      ag[g, ] <- ag1
      a <- m2
    }
  }
  ag
}


# ---------------------------------------------------------------------------
# RANMUL_r
# ---------------------------------------------------------------------------

#' Compute random-mating multipliers and allele copy counts
#'
#' Translates the RANMUL subroutine from de Silva et al. (2005).  For each
#' genotype in \code{ag}, computes:
#' \itemize{
#'   \item \code{rmul[g]}: the multinomial coefficient
#'     \eqn{m2! / \prod_a (\text{count of allele } a \text{ in genotype } g)!},
#'     which is the multiplier for the genotype frequency under random mating.
#'   \item \code{arep[g, a]}: the number of copies of allele \eqn{a} in
#'     genotype \eqn{g}.
#' }
#'
#' @param ng Integer; number of genotypes (rows of \code{ag}).
#' @param na1 Integer; number of alleles including the null allele.
#' @param ag Integer matrix of dimensions \code{ng} x \code{m2} as returned
#'   by \code{\link{GENLIST_r}}.
#' @param m2 Integer; ploidy level.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{rmul}}{Integer vector of length \code{ng}: multinomial
#'       coefficients.}
#'     \item{\code{arep}}{Integer matrix of dimensions \code{ng} x \code{na1}:
#'       allele copy counts.}
#'   }
#'
#' @seealso \code{\link{GENLIST_r}}, \code{\link{SELFMAT_r}}
#'
#' @export
RANMUL_r <- function(ng, na1, ag, m2) {
  # arep: count of each allele in each genotype via tabulate
  arep <- t(apply(ag, 1L, tabulate, nbins = na1))
  storage.mode(arep) <- "integer"

  # rmul: multinomial coefficient = m2! / prod(count_a!)
  rmul <- as.integer(round(exp(
    lfactorial(m2) - rowSums(lfactorial(arep))
  )))

  list(rmul = rmul, arep = arep)
}


# ---------------------------------------------------------------------------
# SELFMAT_r
# ---------------------------------------------------------------------------

#' Build the polyploid selfing matrix
#'
#' Translates the SELFMAT subroutine from de Silva et al. (2005).  Constructs
#' the \code{ng} x \code{ng} integer selfing matrix \code{smat}, where entry
#' \code{smat[g1, g2]} counts the number of gamete-pair combinations from
#' selfing parent \code{g1} that produce offspring genotype \code{g2}.
#'
#' The matrix is divided by \code{(G_r(m - 1, m + 1))^2} (= \code{choose(m2,
#' m)^2}) inside \code{de_silva_freq_EM} to convert counts to proportions
#' (this divisor is precomputed as \code{smatdiv}).
#'
#' @param ng Integer; number of genotypes.
#' @param na1 Integer; number of alleles including the null allele.
#' @param ag Integer matrix (\code{ng} x \code{m2}) of all genotypes, as
#'   returned by \code{\link{GENLIST_r}}.
#' @param m2 Integer; ploidy level (must be even; \code{m = m2 / 2}).
#'
#' @return An integer matrix of dimensions \code{ng} x \code{ng}.
#'
#' @seealso \code{\link{GENLIST_r}}, \code{\link{RANMUL_r}},
#'   \code{\link{INDEXG_r}}
#'
#' @export
SELFMAT_r <- function(ng, na1, ag, m2) {
  m    <- m2 %/% 2L
  smat <- matrix(0L, nrow = ng, ncol = ng)

  # Pre-compute all gamete position combinations (m positions from 1:m2)
  gam_pos <- combn(m2, m)          # m x n_gam matrix
  n_gam   <- ncol(gam_pos)

  # Pre-compute constant for batch INDEXG
  G_top <- as.integer(choose(na1 + m2 - 1L, m2))  # G_r(m2-1, na1)

  # Expand gamete pair indices (all n_gam^2 ordered pairs)
  idx1 <- rep(seq_len(n_gam), each  = n_gam)
  idx2 <- rep(seq_len(n_gam), times = n_gam)

  for (g in seq_len(ng)) {
    parent <- ag[g, ]

    # Extract alleles for each gamete: m x n_gam matrix
    gam_alleles <- matrix(parent[gam_pos], nrow = m, ncol = n_gam)

    # Build all offspring: concatenate gamete pairs, then sort each row
    # off is n_gam^2 x m2
    off <- cbind(
      t(gam_alleles[, idx1, drop = FALSE]),
      t(gam_alleles[, idx2, drop = FALSE])
    )
    off <- t(apply(off, 1L, sort.int, method = "radix"))

    # Batch INDEXG: vectorized over all offspring rows
    x <- 1L + off[, m2] - off[, m2 - 1L]
    if (m2 >= 3L) {
      for (q in 1L:(m2 - 2L)) {
        x <- x +
          as.integer(choose(na1 + 1L - off[, m2 - q - 1L] + q, q + 1L)) -
          as.integer(choose(na1 + 1L - off[, m2 - q]     + q, q + 1L))
      }
    }
    x <- x + G_top -
      as.integer(choose(na1 + 1L - off[, 1L] + m2 - 1L, m2))

    # Accumulate offspring counts into smat row
    smat[g, ] <- tabulate(x, nbins = ng)
  }

  storage.mode(smat) <- "integer"
  smat
}
