# polyfreq-package.R
#
# Package-level roxygen directives: registers the compiled shared library
# produced by src/selfmat_cone.cpp so that .selfmat_cone_cpp() is available
# at runtime.

#' polyfreq: Fast De Silva Allele Frequency Estimation for Polyploids
#'
#' Provides an optimized implementation of the de Silva et al. (2005) EM
#' algorithm for estimating allele frequencies in polyploid populations.
#' Pure-R vectorized code is combined with an optional Rcpp + OpenMP back-end
#' for the computationally dominant sparse selfing matrix step.
#'
#' @docType package
#' @name polyfreq-package
#' @aliases polyfreq
#'
#' @useDynLib polyfreq, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
