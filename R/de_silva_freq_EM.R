# de_silva_freq_EM.R
#
# Main exported function: de_silva_freq_EM
# Pure-R EM-algorithm implementation of De Silva et al. (2005).
#
# Internal closure helpers (not exported):
#   .indexf   — maps a phenotype to its integer index (INDEXP equivalent)
#   .fenlist  — enumerates all phenotypes (PHENLIST equivalent)
#   .convmat  — builds phenotype-to-genotype conversion matrix C
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
#' and produces numerically identical results while replacing all C/C++
#' subroutines (\code{G}, \code{GENLIST}, \code{RANMUL}, \code{SELFMAT}) with
#' vectorized R code.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Precompute per-locus structures: genotype list (\code{ag}),
#'     multinomial multipliers (\code{rmul}), allele copy matrix (\code{arep}),
#'     selfing matrix (\code{smatt}), conversion matrix (\code{cmat}), and
#'     observed phenotype frequencies (\code{pp}).
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
#' @importFrom Matrix Matrix solve
#' @importFrom cli cli_abort
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
                             maxiter  = 1e4) {

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

  # --- Main loop: per locus, per population ----------------------------------

  for (L in loci) {
    cat(paste("Starting", L), sep = "\n")

    for (pop in pops) {
      cat(paste("Starting", L, pop), sep = "\n")

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
      ng <- na1
      for (j in 2L:m2) {
        ng <- ng * (na1 + j - 1L) / j
      }
      ng <- as.integer(ng)

      # Precompute structures
      ag   <- GENLIST_r(ng, na1, m2)
      temp <- .fenlist(na, m2)
      af   <- temp[[1L]]
      naf  <- temp[[2L]]
      nf   <- length(naf)
      temp <- RANMUL_r(ng, na1, ag, m2)
      rmul <- temp$rmul
      arep <- temp$arep
      smatt <- SELFMAT_r(ng, na1, ag, m2) / smatdiv
      cmat  <- .convmat(ng, nf, na1, ag, indexf)

      # Observed phenotype frequencies
      pp <- rep(0, nf)
      for (s in psamples) {
        phenotype <- sort(unique(Genotype(object, s, L)))
        phenotype <- match(phenotype, alleles)
        f <- indexf(length(phenotype), phenotype, na)
        pp[f] <- pp[f] + 1
      }
      pp <- pp / sum(pp)

      # Initial allele frequencies
      p1     <- rep(0, na1)
      p1[na1] <- initNull[L]
      subInitFreq <- subInitFreq * (1 - initNull[L]) / sum(subInitFreq)
      for (a in alleles) {
        p1[match(a, alleles)] <- subInitFreq[1L, paste(L, a, sep = ".")]
      }

      # --- EM algorithm ------------------------------------------------------

      converge <- 0L
      niter    <- 1L
      rmul_d   <- as.numeric(rmul)
      arep_t   <- t(arep)                       # na1 x ng (pre-transpose)

      # Pre-compute (I - s*A)^{-1} — constant across iterations
      s3_inv <- solve(diag(nrow = ng) - self * smatt)

      while (converge == 0L) {

        # -- E-step: expected genotype frequencies ---

        # Allele frequency vector (null = complement)
        pa      <- as.numeric(p1)
        pa[na1] <- 1 - sum(pa[seq_len(na)])

        # Vectorized rvec: rmul[g] * prod(pa[ag[g, ]])
        pa_mat <- matrix(pa[ag], nrow = ng, ncol = m2)
        log_pa <- log(pa_mat)
        log_pa[!is.finite(log_pa)] <- -744  # floor for pa == 0
        rvec <- rmul_d * exp(rowSums(log_pa))

        # Equilibrium genotype probabilities: (1-s) * inv * rvec
        gprob <- (1 - self) * (s3_inv %*% rvec)

        # Eq. (12): distribute genotype probs across phenotypes
        xx1      <- t(t(cmat) * as.numeric(gprob))   # nf x ng
        xx2      <- rowSums(xx1)                      # nf
        xx2_inv  <- ifelse(xx2 > 0, 1 / xx2, 0)
        xx3      <- xx1 * xx2_inv                     # nf x ng (column recycling)
        EP       <- crossprod(xx3, pp)                 # ng x 1

        # -- M-step: update allele frequencies ---
        p2 <- as.numeric(arep_t %*% EP) / m2

        # -- Convergence check ---
        pB <- p1 + p2
        pT <- p1 - p2
        keep <- abs(pB) > 1e-14
        pT <- pT[keep]
        pB <- pB[keep]

        if (length(pB) == 0L || sum(abs(pT) / pB) <= tol) {
          converge <- 1L
        }
        if (niter >= maxiter) {
          converge <- 1L
        }

        niter <- niter + 1L
        p1 <- p2
      }

      # Write final frequencies to output
      for (a in alleles) {
        finalfreq[pop, match(paste(L, a, sep = "."), names(finalfreq))] <-
          p2[match(a, alleles)]
      }
      finalfreq[pop, match(paste(L, "null", sep = "."), names(finalfreq))] <-
        p2[na1]

      cat(paste(niter - 1L, "repetitions for", L, pop), sep = "\n")
    }
  }

  finalfreq
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
