library(polyfreq)
library(polysat)
data_path <- system.file("testdata", "simgen.Rdata", package = "polyfreq")
load(data_path)

# Force sparse path on simgen (small dataset) to verify correctness fix
cat("Running dense path...\n")
res_dense  <- de_silva_freq_EM(simgen, self = 0.1, method = "dense")
cat("Running sparse path...\n")
res_sparse <- de_silva_freq_EM(simgen, self = 0.1, method = "sparse")

# Check freq sums for each locus
freq_loci <- Loci(simgen)
cat("\nFrequency sums by locus (should all be 1.000000):\n")
all_ok <- TRUE
for (L in freq_loci) {
  cols_d <- grep(paste0("^", L, "\\."), names(res_dense))
  cols_s <- grep(paste0("^", L, "\\."), names(res_sparse))
  ds <- rowSums(res_dense[, cols_d, drop=FALSE])
  ss <- rowSums(res_sparse[, cols_s, drop=FALSE])
  ok_d <- all(abs(ds - 1) < 1e-6)
  ok_s <- all(abs(ss - 1) < 1e-6)
  cat(sprintf("  %-8s  dense=[%.6f, %.6f] OK=%s  sparse=[%.6f, %.6f] OK=%s\n",
              L, ds[1], ds[2], ok_d, ss[1], ss[2], ok_s))
  if (!ok_d || !ok_s) all_ok <- FALSE
}

# Max diff between dense and sparse
common_cols <- intersect(names(res_dense), names(res_sparse))
max_diff <- max(abs(res_dense[, common_cols] - res_sparse[, common_cols]))
cat(sprintf("\nMax abs diff dense vs sparse: %.2e\n", max_diff))
cat(sprintf("All freq sums OK: %s\n", all_ok))
