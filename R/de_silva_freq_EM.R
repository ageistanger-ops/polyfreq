# de_silva_freq_EM.R
#
# Main exported function: de_silva_freq_EM
# Pure-R EM-algorithm implementation of De Silva et al. (2005).
#
# Internal helpers (not exported):
#   .make_indexf   — maps a phenotype to its integer index (INDEXP equivalent)
#   .fenlist       — enumerates all phenotypes (PHENLIST equivalent)
#   .convmat       — builds phenotype-to-genotype conversion matrix C (dense path)
#   .build_cone    — builds union phenotype cone (sparse path)
#   .neumann_gprob — Neumann series solver for (I - s*smatt)^{-1} * rvec
#   .one_locus_pop — computes EM for a single (locus, population) pair
#
# Reference: De Silva HN, Hall AJ, Rikkerink E, McNeilage MA, Fraser LG (2005)
#   Estimation of allele frequencies in polyploids under certain patterns of
#   inheritance.  Heredity 95:327-334.  doi:10.1038/sj.hdy.6800728


# ---------------------------------------------------------------------------
# de_silva_freq_EM
# ---------------------------------------------------------------------------

#' EM estimation of allele frequencies in polyploids (de Silva method)
#'
#' An optimized pure-R implementation of the EM algorithm described in
#' De Silva et al. (2005) for estimating allele frequencies in polyploid
#' populations under polysomic inheritance with a mixed mating system
#' (selfing rate \code{self}, random mating rate \code{1 - self}).
#'
#' This function is a drop-in replacement for \code{\link[polysat]{deSilvaFreq}}
#' and produces numerically identical results (dense path) while replacing all
#' C/C++ subroutines (\code{G}, \code{GENLIST}, \code{RANMUL}, \code{SELFMAT})
#' with vectorized R code.
#'
#' For large genotype spaces (\code{ng > 50,000}), a sparse code path is
#' activated automatically.  It restricts computation to the \emph{union
#' phenotype cone} — the set of all genotypes whose non-null alleles are a
#' subset of the alleles observed in at least one sample — and solves the
#' selfing equilibrium via a Neumann series instead of a dense matrix inverse.
#' Because every allele outside the cone has zero initial frequency, this
#' gives the same result as the full-\code{ng} computation.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Precompute per-locus structures: genotype list (\code{ag}),
#'     multinomial multipliers (\code{rmul}), allele copy matrix (\code{arep}),
#'     selfing matrix (\code{smatt}), and observed phenotype frequencies
#'     (\code{pp}).
#'   \item \strong{E-step}: Compute expected equilibrium genotype frequencies
#'     via \eqn{E[P] = (1-s)(I - s A^\top)^{-1} R} (Eq. 6), then distribute
#'     these across phenotypes to obtain expected genotype counts \code{EP}.
#'   \item \strong{M-step}: Update allele frequencies as
#'     \eqn{p_2 = A^\top EP / m2}.
#'   \item Repeat until convergence:
#'     \eqn{\sum_i |p_1^{(i)} - p_2^{(i)}| / (p_1^{(i)} + p_2^{(i)}) \le tol}.
#' }
#'
#' @param object A \code{\link[polysat]{genambig}} (or \code{genbinary})
#'   object containing the polyploid genotype data.
#' @param self Numeric scalar in \eqn{[0, 1]}; the selfing rate (required).
#' @param samples Character vector of sample names to include.  Defaults to
#'   all samples in \code{object}.
#' @param loci Character vector of locus names to include.  Defaults to all
#'   loci in \code{object}.
#' @param initNull Numeric scalar or named vector; initial null-allele
#'   frequency.  A single value is recycled across all loci.  Default is
#'   \code{0.15}.
#' @param initFreq Data frame of initial allele frequencies as returned by
#'   \code{\link[polysat]{simpleFreq}}.  Must contain a \code{"Genomes"}
#'   column.  Defaults to \code{polysat::simpleFreq(object[samples, loci])}.
#' @param tol Numeric; convergence tolerance (default \code{1e-8}).
#' @param maxiter Integer; maximum number of EM iterations (default
#'   \code{10000}).
#' @param method Character; code path selector: \code{"auto"} (default)
#'   switches to the sparse path when \code{ng > 50000}, \code{"dense"}
#'   always uses the dense path (may fail for large \code{ng}), \code{"sparse"}
#'   always uses the sparse path.
#' @param nthreads Integer; number of parallel workers for the locus x
#'   population loop (default \code{1L} = sequential).  Values \code{> 1}
#'   require the \pkg{future} and \pkg{future.apply} packages.  When
#'   \code{nthreads > 1}, the Rcpp sparse selfmat uses 1 OpenMP thread per
#'   worker to avoid CPU over-subscription; when \code{nthreads == 1}, the
#'   Rcpp selfmat uses all available cores via OpenMP.
#'
#' @return A data frame with one row per population and one column per allele
#'   (named \code{"Locus.allele"}) plus one null column per locus (named
#'   \code{"Locus.null"}).  A \code{"Genomes"} column (genome count) is
#'   included as the first column, matching the format of
#'   \code{\link[polysat]{deSilvaFreq}}.
#'
#' @seealso \code{\link[polysat]{deSilvaFreq}}, \code{\link{GENLIST_r}},
#'   \code{\link{RANMUL_r}}, \code{\link{SELFMAT_r}}, \code{\link{G_r}}
#'
#' @references
#'   De Silva HN, Hall AJ, Rikkerink E, McNeilage MA, Fraser LG (2005).
#'   Estimation of allele frequencies in polyploids under certain patterns of
#'   inheritance. \emph{Heredity} 95:327-334.
#'   \doi{10.1038/sj.hdy.6800728}
#'
#' @importFrom polysat Samples Loci Ploidies PopInfo PopNames Genotype
#'   isMissing simpleFreq genbinary.to.genambig
#' @importFrom Matrix Matrix solve sparseMatrix
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#'
#' @examples
#' \dontrun{
#' library(polysat)
#' data_path <- system.file("testdata", "simgen.Rdata", package = "polyfreq")
#' load(data_path)
#' result <- de_silva_freq_EM(simgen, self = 0.1)
#' head(result)
#' }
#'
#' @export
de_silva_freq_EM <- function(object, self,
                             samples  = Samples(object),
                             loci     = Loci(object),
                             initNull = 0.15,
                             initFreq = simpleFreq(object[samples, loci]),
                             tol      = 1e-8,
                             maxiter  = 1e4,
                             method   = c("auto", "dense", "sparse"),
                             nthreads = 1L) {

  method   <- match.arg(method)
  nthreads <- as.integer(nthreads)

  # --- Input validation (matching polysat::deSilvaFreq) ----------------------

  if (missing(self))
    cli_abort("Selfing rate required.")

  if (!"Genomes" %in% names(initFreq))
    cli_abort("initFreq must have single {.field Genomes} column.")

  # Convert genbinary to genambig if necessary
  if ("genbinary" %in% class(object))
    object <- genbinary.to.genambig(object, samples, loci)

  # Get and validate ploidy
  m2 <- unique(as.vector(Ploidies(object, samples, loci)))
  if (length(m2) != 1L)
    cli_abort("Only one ploidy allowed. Try running subsets of data one ploidy at a time.")
  if (is.na(m2))
    cli_abort("Function requires information in Ploidies slot.")
  if (m2 %% 2L != 0L)
    cli_abort("Ploidy must be even.")

  # Check and set up initNull
  if (!length(initNull) %in% c(1L, length(loci)))
    cli_abort("Need single value for initNull or one value per locus.")
  if (length(initNull) == 1L)
    initNull <- rep(initNull, times = length(loci))
  if (is.null(names(initNull)))
    names(initNull) <- loci

  m <- m2 %/% 2L

  # Get populations
  pops <- PopNames(object)[unique(PopInfo(object)[samples])]
  if (!identical(pops, row.names(initFreq)))
    cli_abort("Population names in initFreq don't match those in object.")

  # --- Output data frame setup (matching polysat format exactly) -------------

  finalfreq <- data.frame(
    row.names = pops,
    Genomes   = initFreq$Genomes,
    matrix(0, nrow = length(pops),
           ncol = dim(initFreq)[2] - 1L + length(loci),
           dimnames = list(NULL, sort(c(
             paste(loci, ".null", sep = ""),
             names(initFreq)[2:dim(initFreq)[2]]
           ))))
  )

  # Divisor for the selfing matrices
  smatdiv <- (G_r(m - 1L, m + 1L))^2

  # Create indexf closure
  indexf <- .make_indexf(m2)

  # --- Build task list: one entry per (locus, population) pair ---------------

  locus_pop_df <- expand.grid(L = loci, pop = pops, stringsAsFactors = FALSE)
  tasks        <- split(locus_pop_df, seq_len(nrow(locus_pop_df)))

  # Shared read-only arguments passed to every worker
  shared <- list(
    object   = object,
    self     = self,
    initFreq = initFreq,
    initNull = initNull,
    tol      = tol,
    maxiter  = maxiter,
    method   = method,
    m2       = m2,
    smatdiv  = smatdiv,
    samples  = samples
  )

  # Number of OpenMP threads for the Rcpp selfmat within each worker.
  # Partition cores proportionally so that (nthreads R workers) ×
  # (cpp_threads OpenMP threads) ≈ total physical cores.
  # nthreads=1  → all cores to OpenMP (e.g. 12 threads on a 12-core machine)
  # nthreads=4  → 3 OpenMP threads per worker  (12 / 4)
  # nthreads=12 → 1 OpenMP thread per worker   (12 / 12)
  cpp_threads <- max(1L, parallel::detectCores() %/% nthreads)

  # --- Dispatch ---------------------------------------------------------------

  worker_fn <- function(pair) {
    do.call(polyfreq:::.one_locus_pop,
            c(pair,
              shared,
              list(quiet       = nthreads > 1L,
                   cpp_threads = cpp_threads)))
  }

  if (nthreads == 1L) {
    # Sequential — identical to old double-loop behaviour
    results <- lapply(tasks, function(pair) {
      do.call(.one_locus_pop,
              c(pair,
                shared,
                list(quiet       = FALSE,
                     cpp_threads = cpp_threads)))
    })
  } else {
    future::plan(future::multisession, workers = nthreads)
    on.exit(future::plan(future::sequential), add = TRUE)
    results <- future.apply::future_lapply(
      tasks,
      worker_fn,
      future.packages = "polyfreq",
      future.seed     = TRUE
    )
  }

  # --- Assemble results into finalfreq ----------------------------------------

  for (res in results) {
    pop     <- res$pop
    L       <- res$L
    p2      <- res$p2
    alleles <- res$alleles
    na1     <- res$na1

    for (a in alleles) {
      finalfreq[pop, match(paste(L, a, sep = "."), names(finalfreq))] <-
        p2[match(a, alleles)]
    }
    finalfreq[pop, match(paste(L, "null", sep = "."), names(finalfreq))] <-
      p2[na1]
  }

  finalfreq
}


# ---------------------------------------------------------------------------
# .one_locus_pop
# ---------------------------------------------------------------------------
# Internal worker: runs the EM algorithm for a single (locus, population)
# pair.  Called sequentially or in parallel via future_lapply.
#
# Arguments:
#   L, pop      character   locus name and population name
#   object      genambig    the full dataset (read-only)
#   self        numeric     selfing rate
#   initFreq    data.frame  initial allele frequencies (from simpleFreq)
#   initNull    named num   initial null frequency per locus
#   tol         numeric     convergence tolerance
#   maxiter     integer     max EM iterations
#   method      character   "auto" | "dense" | "sparse"
#   m2          integer     ploidy
#   smatdiv     integer     selfing matrix divisor (choose(m2,m)^2)
#   samples     character   sample names in scope
#   quiet       logical     suppress cli progress messages (TRUE in workers)
#   cpp_threads integer     OpenMP threads for Rcpp selfmat
#
# Returns a named list: list(L, pop, p2, alleles, na1)
.one_locus_pop <- function(L, pop, object, self, initFreq, initNull,
                            tol, maxiter, method, m2, smatdiv, samples,
                            quiet = FALSE, cpp_threads = 1L) {

  m      <- m2 %/% 2L
  indexf <- .make_indexf(m2)

  # Samples in this population with non-missing data at this locus
  psamples <- Samples(object, populations = pop)[
    !isMissing(object, Samples(object, populations = pop), L)
  ]
  psamples <- psamples[psamples %in% samples]

  # Extract initial allele frequencies for this locus-population
  subInitFreq <- initFreq[pop, grep(paste("^", L, "\\.", sep = ""),
                                    names(initFreq))]

  # Identify alleles with non-zero initial frequency
  templist <- names(subInitFreq)[subInitFreq != 0]
  templist <- strsplit(templist, split = ".", fixed = TRUE)
  alleles  <- vapply(templist, function(x) as.integer(x[2L]), integer(1L))
  alleles  <- sort(alleles)
  na  <- length(alleles)
  na1 <- na + 1L

  # Number of genotypes: choose(na1 + m2 - 1, m2)
  ng <- as.numeric(choose(na1 + m2 - 1L, m2))

  # Decide which code path to use
  use_sparse <- switch(method,
    "dense"  = FALSE,
    "sparse" = TRUE,
    "auto"   = (ng > 50000)
  )

  # Initial allele frequencies
  p1      <- rep(0, na1)
  p1[na1] <- initNull[L]
  subInitFreq <- subInitFreq * (1 - initNull[L]) / sum(subInitFreq)
  for (a in alleles) {
    p1[match(a, alleles)] <- subInitFreq[1L, paste(L, a, sep = ".")]
  }

  # --- Progress: start -------------------------------------------------------
  if (!quiet) {
    path_label <- if (use_sparse) "sparse" else "dense "
    cli_alert_info(
      "[{L} / {pop}]  {path_label} path  (ng = {format(ng, big.mark = ',', scientific = FALSE)})"
    )
  }

  if (!use_sparse) {
    # ========================================================================
    # DENSE PATH — bit-identical to polysat::deSilvaFreq
    # ========================================================================

    ng_int <- as.integer(ng)

    # Precompute structures
    ag   <- GENLIST_r(ng_int, na1, m2)
    temp <- .fenlist(na, m2)
    af   <- temp[[1L]]
    naf  <- temp[[2L]]
    nf   <- length(naf)
    temp <- RANMUL_r(ng_int, na1, ag, m2)
    rmul <- temp$rmul
    arep <- temp$arep
    smatt <- SELFMAT_r(ng_int, na1, ag, m2) / smatdiv
    cmat  <- .convmat(ng_int, nf, na1, ag, indexf)

    # Observed phenotype frequencies
    pp <- rep(0, nf)
    for (s in psamples) {
      phenotype <- sort(unique(Genotype(object, s, L)))
      phenotype <- match(phenotype, alleles)
      f <- indexf(length(phenotype), phenotype, na)
      pp[f] <- pp[f] + 1
    }
    pp <- pp / sum(pp)

    # --- EM algorithm --------------------------------------------------

    converge <- 0L
    niter    <- 1L
    rmul_d   <- as.numeric(rmul)
    arep_t   <- t(arep)                     # na1 x ng (pre-transpose)

    # Pre-compute (I - s*A)^{-1} — constant across iterations
    s3_inv <- solve(diag(nrow = ng_int) - self * smatt)

    while (converge == 0L) {

      # E-step: expected genotype frequencies
      pa      <- as.numeric(p1)
      pa[na1] <- 1 - sum(pa[seq_len(na)])

      pa_mat <- matrix(pa[ag], nrow = ng_int, ncol = m2)
      log_pa <- log(pa_mat)
      log_pa[!is.finite(log_pa)] <- -744
      rvec <- rmul_d * exp(rowSums(log_pa))

      gprob <- (1 - self) * (s3_inv %*% rvec)

      # Eq. (12): distribute genotype probs across phenotypes
      xx1      <- t(t(cmat) * as.numeric(gprob))
      xx2      <- rowSums(xx1)
      xx2_inv  <- ifelse(xx2 > 0, 1 / xx2, 0)
      xx3      <- xx1 * xx2_inv
      EP       <- crossprod(xx3, pp)

      # M-step: update allele frequencies
      p2 <- as.numeric(arep_t %*% EP) / m2

      # Convergence check
      pB   <- p1 + p2
      pT   <- p1 - p2
      keep <- abs(pB) > 1e-14
      pT   <- pT[keep]
      pB   <- pB[keep]

      if (length(pB) == 0L || sum(abs(pT) / pB) <= tol) {
        converge <- 1L
      }
      if (niter >= maxiter) {
        converge <- 1L
      }

      niter <- niter + 1L
      p1 <- p2
    }

  } else {
    # ========================================================================
    # SPARSE PATH — union phenotype cone + Neumann series
    # ========================================================================

    # --- Build phenotype cone ------------------------------------------
    cone        <- .build_cone(psamples, object, L, alleles, na, na1,
                               m2, indexf)
    ag_c        <- cone$ag_cone
    n_c         <- cone$n_cone
    cph         <- cone$compact_ph_idx   # compact phenotype index (1..n_uph)
    unique_ph   <- cone$unique_ph        # original indexf value per compact group
    n_uph       <- cone$n_uph            # number of unique phenotypes in cone

    if (!quiet)
      cli_alert_info(
        "  cone: {format(n_c, big.mark = ',')} genotypes, {n_uph} phenotypes"
      )

    # Precompute RANMUL on cone genotypes
    temp   <- RANMUL_r(n_c, na1, ag_c, m2)
    rmul_c <- as.numeric(temp$rmul)
    arep_c <- temp$arep

    # Build sparse selfing matrix on cone (one-time cost)
    smatt_c <- .selfmat_cone_batch(ag_c, na1, m2,
                                   nthreads = cpp_threads) / smatdiv

    # Observed phenotype frequencies — compact vector of length n_uph.
    ph_to_compact <- function(f) match(f, unique_ph)

    pp_compact <- numeric(n_uph)
    for (s in psamples) {
      phenotype <- sort(unique(Genotype(object, s, L)))
      phenotype <- match(phenotype, alleles)
      f         <- indexf(length(phenotype), phenotype, na)
      j         <- ph_to_compact(f)
      if (!is.na(j)) pp_compact[j] <- pp_compact[j] + 1L
    }
    pp_compact <- pp_compact / sum(pp_compact)

    # --- EM algorithm (sparse) -----------------------------------------

    converge  <- 0L
    niter     <- 1L
    arep_c_t  <- t(arep_c)          # na1 x n_c (pre-transpose)

    while (converge == 0L) {

      # E-step: random-mating genotype frequencies
      pa      <- as.numeric(p1)
      pa[na1] <- 1 - sum(pa[seq_len(na)])

      pa_mat <- matrix(pa[ag_c], nrow = n_c, ncol = m2)
      log_pa <- log(pa_mat)
      log_pa[!is.finite(log_pa)] <- -744
      rvec <- rmul_c * exp(rowSums(log_pa))

      # Selfing equilibrium via Neumann series
      gprob <- .neumann_gprob(smatt_c, rvec, self)

      # Distribute genotype probs across compact phenotype groups.
      xx2_compact <- as.numeric(
        rowsum(matrix(gprob, ncol = 1L), cph, reorder = TRUE)
      )
      xx2 <- xx2_compact[cph]    # expand back to length n_c

      # Expected genotype counts
      EP <- gprob * pp_compact[cph] / xx2
      EP[!is.finite(EP)] <- 0

      # M-step: update allele frequencies
      p2 <- as.numeric(arep_c_t %*% EP) / m2

      # Convergence check
      pB   <- p1 + p2
      pT   <- p1 - p2
      keep <- abs(pB) > 1e-14
      pT   <- pT[keep]
      pB   <- pB[keep]

      if (length(pB) == 0L || sum(abs(pT) / pB) <= tol) {
        converge <- 1L
      }
      if (niter >= maxiter) {
        converge <- 1L
      }

      niter <- niter + 1L
      p1 <- p2
    }

  }

  # --- Progress: done --------------------------------------------------------
  if (!quiet)
    cli_alert_success("[{L} / {pop}]  {niter - 1L} iterations")

  list(L = L, pop = pop, p2 = p2, alleles = alleles, na1 = na1)
}


# ---------------------------------------------------------------------------
# Internal closure helpers (not exported)
# ---------------------------------------------------------------------------

# .make_indexf --------------------------------------------------------------
# Maps a phenotype (sorted allele vector) to its integer index.
# Pure-R translation of the INDEXF subroutine from de Silva et al. (2005).
#
# Arguments:
#   m    integer  number of distinct alleles in this phenotype
#   af1  integer vector  sorted allele indices (length >= m, only first m used)
#   na   integer  number of non-null alleles at this locus
#
# Note: this is defined as a closure factory so that m2 is captured from
# the enclosing scope, matching the polysat design where indexf is defined
# inside deSilvaFreq.
.make_indexf <- function(m2) {
  function(m, af1, na) {
    x <- 1L
    if (m == 0L) return(x)
    if (m == 1L) {
      x <- x + af1[1L]
    } else if (m > 1L) {
      for (q in 1L:(m - 1L)) {
        x <- x + G_r(q - 1L, na - q + 1L)
      }
      x <- x + G_r(m - 1L, na + 1L - m) - G_r(m - 1L, na + 2L - m - af1[1L])
      if (m > 2L) {
        for (q in 1L:(m - 2L)) {
          x <- x + G_r(q, na - q - af1[m - q - 1L]) -
                   G_r(q, na + 1L - q - af1[m - q])
        }
      }
      x <- x + af1[m] - af1[m - 1L]
    }
    as.integer(x)
  }
}


# .fenlist -----------------------------------------------------------------
# Enumerates all phenotypes (observable allele subsets, ignoring copy number).
# Pure-R translation of the FENLIST subroutine from de Silva et al. (2005).
#
# Arguments:
#   na   integer  number of non-null alleles at this locus
#   m2   integer  ploidy level
#
# Returns a list(af, naf) (unnamed, positional) where:
#   af   integer matrix (nf x m2) — each row is one phenotype (padded with 0)
#   naf  integer vector (nf)      — number of distinct alleles per phenotype
.fenlist <- function(na, m2) {
  # First phenotype is null (all zeros)
  af1 <- rep(0L, m2)
  af  <- matrix(af1, nrow = 1L)
  naf <- 0L

  for (m_size in seq_len(min(m2, na))) {
    # Initial combination of m_size alleles from 1:na
    af1[seq_len(m_size)] <- seq_len(m_size)
    if (m_size < m2) af1[(m_size + 1L):m2] <- 0L

    naf <- c(naf, m_size)
    af  <- rbind(af, af1)

    # Iterate remaining combinations via backtracking
    a <- m_size
    while (a > 0L) {
      if (af1[a] == (na + a - m_size)) {
        a <- a - 1L
      } else {
        af1[a] <- af1[a] + 1L
        if (a < m_size) {
          for (a1 in (a + 1L):m_size) af1[a1] <- af1[a1 - 1L] + 1L
        }
        naf <- c(naf, m_size)
        af  <- rbind(af, af1)
        a <- m_size
      }
    }
  }

  storage.mode(af) <- "integer"
  list(af, as.integer(naf))
}


# .convmat -----------------------------------------------------------------
# Builds the nf x ng phenotype-to-genotype conversion matrix C.
# Pure-R translation of the CONVMAT subroutine.
#
# Arguments:
#   ng      integer                number of genotypes
#   nf      integer                number of phenotypes
#   na1     integer                number of alleles including null
#   ag      integer matrix (ng x m2) genotype list from GENLIST_r
#   indexf  function               closure returned by .make_indexf(m2)
#
# Returns an integer matrix (nf x ng); cmat[f, g] = 1 iff genotype g
# produces phenotype f.
.convmat <- function(ng, nf, na1, ag, indexf) {
  na   <- na1 - 1L
  m2   <- ncol(ag)
  af1  <- rep(0L, m2)
  cmat <- matrix(0L, nrow = nf, ncol = ng)

  for (g in seq_len(ng)) {
    ag1 <- ag[g, ]

    if (ag1[1L] == na1) {
      # Homozygous null genotype
      naf1 <- 0L
    } else {
      naf1      <- 1L
      af1[naf1] <- ag1[1L]
      for (a in 2L:m2) {
        if (ag1[a] == na1) break
        if (ag1[a] > af1[naf1]) {
          naf1      <- naf1 + 1L
          af1[naf1] <- ag1[a]
        }
      }
    }

    # Zero out unused positions
    if (naf1 < m2) {
      af1[(naf1 + 1L):m2] <- 0L
    }

    f <- indexf(naf1, af1, na)
    cmat[f, g] <- 1L
  }

  cmat
}


# .build_cone --------------------------------------------------------------
# Builds the union phenotype cone for the sparse code path.
#
# The cone = union over all observed phenotypes P of:
#   { all genotypes whose non-null alleles form any subset of P }
# = union over all P of
#   { all m2-multisets from alleles(P) ∪ {null} }
#
# This set is closed under selfing (offspring alleles ⊆ parent alleles ⊆ cone
# alleles), so running the EM within the cone gives the same result as the
# full-ng computation (genotypes outside the cone contribute zero to rvec
# because every allele absent from the cone has simpleFreq = 0).
#
# NOTE: nf (the total number of phenotypes) can be astronomically large for
# octoploids with many alleles (e.g., C(27,8) ≈ 3.5 M for na=27, m2=8).
# We therefore return COMPACT phenotype indices (1..n_uph) based only on the
# phenotype indices actually present in the cone, avoiding any large nf-sized
# vector allocation.
#
# Arguments:
#   psamples  character vector   sample IDs in this population
#   object    genambig           the full dataset
#   L         character          locus name
#   alleles   integer vector     global allele values (length na)
#   na        integer            number of non-null alleles
#   na1       integer            na + 1 (null allele index in local numbering)
#   m2        integer            ploidy
#   indexf    function           closure from .make_indexf(m2)
#
# Returns a named list:
#   ag_cone         integer matrix (n_cone x m2) — cone genotypes
#   n_cone          integer                       — number of cone genotypes
#   compact_ph_idx  integer vector (n_cone)       — compact phenotype index (1..n_uph)
#   unique_ph       integer vector (n_uph)        — original indexf values per compact group
#   n_uph           integer                       — number of unique phenotypes in cone
#   cone_hashes     numeric vector (n_cone)       — polynomial hash per row
#   hash_weights    numeric vector (m2)           — weights used for hashing
.build_cone <- function(psamples, object, L, alleles, na, na1, m2, indexf) {

  # --- Collect unique observed phenotypes -----------------------------------
  # Each entry of unique_phenos is a sorted integer vector of positions in
  # alleles[] (values 1..na, no null).
  # polysat::Genotype() uses S4 dispatch and CANNOT be called from C++,
  # so this loop stays in R.
  pheno_seen <- list()
  for (s in psamples) {
    raw       <- sort(unique(Genotype(object, s, L)))
    local_idx <- match(raw, alleles)      # positions 1..na
    key       <- paste(local_idx, collapse = "-")
    if (!key %in% names(pheno_seen)) {
      pheno_seen[[key]] <- local_idx
    }
  }
  unique_phenos <- unname(pheno_seen)

  # --- C++: GENLIST per phenotype + dedup + INDEXF + compact mapping --------
  # .build_cone_cpp() (src/build_cone.cpp) handles:
  #   1. GENLIST enumeration of m2-multisets for each phenotype closure
  #   2. Global deduplication by polynomial hash
  #   3. INDEXF computation for every unique cone genotype
  #   4. Compact phenotype mapping (sort/unique/match)
  result <- .build_cone_cpp(unique_phenos, na, na1, m2)

  ag_cone <- result$ag_cone
  n_cone  <- nrow(ag_cone)

  # Recompute hash weights and cone hashes with base=na1 (NOT na1+1) to
  # match the convention used by .selfmat_cone_batch() / .selfmat_cone_cpp().
  # These are not consumed by .one_locus_pop() itself but are kept in the
  # return value for interface consistency.
  hash_weights <- as.numeric(na1)^(seq(0L, m2 - 1L))
  cone_hashes  <- as.numeric(ag_cone %*% hash_weights)

  list(
    ag_cone        = ag_cone,
    n_cone         = n_cone,
    compact_ph_idx = result$compact_ph_idx,
    unique_ph      = result$unique_ph,
    n_uph          = length(result$unique_ph),
    cone_hashes    = cone_hashes,
    hash_weights   = hash_weights
  )
}


# .neumann_gprob -----------------------------------------------------------
# Computes gprob = (1-s) * (I - s * smatt)^{-1} * rvec via Neumann series:
#   gprob = (1-s) * sum_{k=0}^{inf} s^k * smatt^k * rvec
#
# Replaces the dense `solve(diag(ng) - self * smatt) %*% rvec`.
# smatt is column-stochastic after dividing by smatdiv, so the series
# converges for any s in [0, 1).  For s = 0, returns rvec unchanged.
#
# Arguments:
#   smatt_sparse  dgCMatrix     sparse selfing matrix (n_cone x n_cone)
#   rvec          numeric       random-mating frequency vector (length n_cone)
#   self          numeric       selfing rate in [0, 1)
#   tol           numeric       stop when s^k * max(|v|) < tol (default 1e-12)
#   maxterms      integer       hard iteration limit (default 200)
#
# Returns a numeric vector of length n_cone.
.neumann_gprob <- function(smatt_sparse, rvec, self,
                           tol      = 1e-12,
                           maxterms = 200L) {
  if (self == 0) return(rvec)

  result <- rvec * (1 - self)
  v      <- rvec
  sk     <- self          # tracks s^k

  for (k in seq_len(maxterms)) {
    v      <- as.numeric(smatt_sparse %*% v)
    result <- result + (1 - self) * sk * v
    sk     <- sk * self
    if (sk * max(abs(v)) < tol) break
  }

  result
}
