library(polyfreq)
library(polysat)

data_path <- system.file("testdata", "hornbeam_genambig.Rdata", package = "polyfreq")
load(data_path)

loci_list <- polysat::Loci(hornbeam_genambig)
cat("Loci:", paste(loci_list, collapse = ", "), "\n")
cat("Ploidy:", unique(as.vector(polysat::Ploidies(hornbeam_genambig))), "\n")
cat("Pops:", paste(polysat::PopNames(hornbeam_genambig), collapse = ", "), "\n\n")

total_start <- proc.time()["elapsed"]

for (L in loci_list) {
  t0 <- proc.time()["elapsed"]
  tryCatch({
    result  <- de_silva_freq_EM(hornbeam_genambig, self = 0.5, loci = L)
    elapsed <- proc.time()["elapsed"] - t0

    freq_pat  <- paste0("^", L, "\\.")
    freq_cols <- grep(freq_pat, names(result))
    freq_sums <- rowSums(result[, freq_cols, drop = FALSE])
    ok        <- all(abs(freq_sums - 1) < 1e-6)
    cat(sprintf("Locus %-8s  %6.1f s  freqs_sum_to_1=%s\n", L, elapsed, ok))
  }, error = function(e) {
    elapsed <- proc.time()["elapsed"] - t0
    cat(sprintf("Locus %-8s  %6.1f s  ERROR: %s\n", L, elapsed,
                conditionMessage(e)))
  })
}

total_elapsed <- proc.time()["elapsed"] - total_start
cat(sprintf("\nTotal: %.1f s\n", total_elapsed))
