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


# ---------------------------------------------------------------------------
# Internal helpers for the sparse code path (not exported)
# ---------------------------------------------------------------------------

# .merge_sorted_pairs -------------------------------------------------------
# Vectorised merge of two sorted m-column integer matrices into a sorted
# 2m-column integer matrix.  Processes N rows simultaneously in 2m passes.
#
# Arguments:
#   g1  integer matrix (N x m)  left sorted halves
#   g2  integer matrix (N x m)  right sorted halves
#   m   integer                 half-ploidy (m2 / 2)
#
# Returns an integer matrix (N x 2m) whose rows are the row-wise merges of
# the sorted halves g1 and g2.
.merge_sorted_pairs <- function(g1, g2, m) {
  n   <- nrow(g1)
  m2  <- 2L * m

  # Append a sentinel column so matrix indexing never goes out of bounds
  sentinel <- rep(.Machine$integer.max, n)
  g1 <- cbind(g1, sentinel)
  g2 <- cbind(g2, sentinel)

  result <- matrix(0L, nrow = n, ncol = m2)
  ia     <- rep(1L, n)
  ib     <- rep(1L, n)
  rows   <- seq_len(n)

  for (k in seq_len(m2)) {
    a_val  <- g1[cbind(rows, ia)]
    b_val  <- g2[cbind(rows, ib)]
    take_a <- a_val <= b_val
    result[, k] <- ifelse(take_a, a_val, b_val)
    ia <- ia + take_a
    ib <- ib + 1L - take_a
  }
  result
}


# .selfmat_cone_batch -------------------------------------------------------
# Builds a sparse selfing matrix for a set of cone genotypes using batch-
# vectorised processing.  Replaces the slow per-parent R loop in SELFMAT_r
# with fully vectorised gamete extraction + merge-sort + hash lookup.
#
# The returned matrix M satisfies M[i, j] = number of gamete-pair combinations
# of parent i that produce offspring j.  Divide by smatdiv in the caller to
# convert to proportions (same convention as SELFMAT_r).
#
# Arguments:
#   ag_cone    integer matrix (n_cone x m2) — cone genotypes
#   na1        integer                       — allele count incl. null
#   m2         integer                       — ploidy (must be even)
#   chunk_size integer                       — parents per batch (default 200)
#
# Returns a Matrix::sparseMatrix of class "dgCMatrix" (n_cone x n_cone).
.selfmat_cone_batch <- function(ag_cone, na1, m2,
                                chunk_size = 200L) {
  n_cone <- nrow(ag_cone)
  m      <- m2 %/% 2L

  # --- Pre-compute once ---------------------------------------------------

  # Gamete position combinations: each column of gam_pos selects m positions
  # from 1:m2 that make up one gamete.  n_gam = C(m2, m) = 70 for octoploid.
  gam_pos <- combn(m2, m)          # m x n_gam
  n_gam   <- ncol(gam_pos)
  # Unordered gamete pairs: sorted_merge(g_a, g_b) = sorted_merge(g_b, g_a),
  # so (i,j) and (j,i) give identical offspring.  Use only i<=j pairs:
  #   - homogametic  (i == i): n_gam pairs,       weight 1
  #   - heterogametic (i < j): C(n_gam,2) pairs,  weight 2
  # Reduces from n_gam^2 = 4900 to C(n_gam+1, 2) = 2485 pairs for octoploid.
  het_pairs <- combn(n_gam, 2L)                  # 2 x C(n_gam, 2)
  idx1    <- c(seq_len(n_gam), het_pairs[1L, ])  # n_pairs total
  idx2    <- c(seq_len(n_gam), het_pairs[2L, ])
  pair_wt <- c(rep(1L, n_gam), rep(2L, ncol(het_pairs)))
  n_pairs <- length(idx1)                         # C(n_gam+1, 2) = 2485

  # Polynomial hash weights: must be > na1 to avoid collisions
  # Use double arithmetic to avoid integer overflow for large na1
  hash_weights <- as.numeric(na1)^(seq(0L, m2 - 1L))

  # Hash for every cone genotype row
  cone_hashes <- as.numeric(ag_cone %*% hash_weights)

  # --- Accumulate sparse triplets -----------------------------------------
  # Collect per-chunk into lists (avoids O(total_nnz^2) from growing vectors).
  n_chunks  <- ceiling(n_cone / chunk_size)
  si_chunks <- vector("list", n_chunks)
  sj_chunks <- vector("list", n_chunks)
  sx_chunks <- vector("list", n_chunks)
  chunk_k   <- 1L

  for (chunk_start in seq(1L, n_cone, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1L, n_cone)
    B         <- chunk_end - chunk_start + 1L   # actual chunk size

    batch <- ag_cone[chunk_start:chunk_end, , drop = FALSE]   # B x m2

    # ---- Step 1: build gamete matrix (B*n_gam) x m ---
    # For each of the B parents, compute all n_gam gametes; stack vertically.
    parent_idx <- rep(seq_len(B), each  = n_gam)   # which parent (1..B)
    gam_col_j  <- rep(seq_len(n_gam), times = B)   # which gamete (1..n_gam)

    # gamete_mat[r, k] = allele at position gam_pos[k, gam_col_j[r]] of parent r
    gamete_mat <- matrix(0L, nrow = B * n_gam, ncol = m)
    for (k in seq_len(m)) {
      gamete_mat[, k] <- batch[cbind(parent_idx,
                                     gam_pos[k, gam_col_j])]
    }

    # ---- Step 2: build offspring pairs (B*n_pairs) x m each ---
    # For each parent b and each ordered pair (idx1[p], idx2[p]):
    #   gamete rows for parent b live at (b-1)*n_gam + 1:n_gam
    base_offset <- rep((seq_len(B) - 1L) * n_gam, each = n_pairs)
    g1_rows     <- base_offset + rep(idx1, times = B)
    g2_rows     <- base_offset + rep(idx2, times = B)

    g1 <- gamete_mat[g1_rows, , drop = FALSE]   # (B*n_pairs) x m
    g2 <- gamete_mat[g2_rows, , drop = FALSE]   # (B*n_pairs) x m

    # ---- Step 3: merge-sort each offspring row ---
    off <- .merge_sorted_pairs(g1, g2, m)        # (B*n_pairs) x m2

    # ---- Step 4: hash offspring and look up cone membership ---
    off_hash <- as.numeric(off %*% hash_weights)
    j_local  <- match(off_hash, cone_hashes)     # NA if not in cone

    # Pair weights for this chunk (1 = homogametic, 2 = heterogametic)
    wt_local <- rep(pair_wt, times = B)

    # ---- Step 5 (fully vectorised): per-parent weighted count → sparse triplets ---
    # Replace the B-iteration inner loop with a single order() + rle() over
    # all valid offspring in the chunk.  This eliminates ~200 R-level function
    # calls per chunk (48 K total) and removes within-chunk c() growth.
    #
    # Encoding: pair_key = global_par * (n_cone + 1) + j_valid.
    # Values fit exactly in double (max ≈ 48400 * 48401 < 2.34e9 < 2^53). ✓
    # Sorted order groups by global_par first, then by j_valid within parent,
    # so rle() on pair_key gives (parent, offspring) pairs with counts.
    valid_mask <- !is.na(j_local)
    if (any(valid_mask)) {
      j_valid    <- j_local[valid_mask]
      wt_valid   <- wt_local[valid_mask]
      par_local  <- rep(seq_len(B), each = n_pairs)[valid_mask]   # 1..B
      global_par <- chunk_start - 1L + par_local

      ord      <- order(global_par, j_valid)
      par_ord  <- global_par[ord]
      j_ord    <- j_valid[ord]
      wt_ord   <- as.numeric(wt_valid[ord])
      pair_key <- as.numeric(par_ord) * as.numeric(n_cone + 1L) + j_ord
      rl       <- rle(pair_key)

      # Sum weights per (parent, offspring) run using cumulative sum trick
      cum_wt   <- cumsum(wt_ord)
      sx_c     <- diff(c(0.0, cum_wt[cumsum(rl$lengths)]))

      si_c <- as.integer(rl$values %/% as.numeric(n_cone + 1L))
      sj_c <- as.integer(rl$values %%  as.numeric(n_cone + 1L))
    } else {
      si_c <- integer(0)
      sj_c <- integer(0)
      sx_c <- numeric(0)
    }

    si_chunks[[chunk_k]] <- si_c
    sj_chunks[[chunk_k]] <- sj_c
    sx_chunks[[chunk_k]] <- sx_c
    chunk_k <- chunk_k + 1L
  }

  si <- unlist(si_chunks, use.names = FALSE)
  sj <- unlist(sj_chunks, use.names = FALSE)
  sx <- unlist(sx_chunks, use.names = FALSE)

  Matrix::sparseMatrix(i    = si,
                       j    = sj,
                       x    = sx,
                       dims = c(n_cone, n_cone))
}
