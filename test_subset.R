library(polyfreq)
library(polysat)

data_path <- system.file("testdata", "hornbeam_genambig.Rdata", package = "polyfreq")
load(data_path)

# --- Identify first population and locus with fewest alleles ----------------
pop1      <- PopNames(hornbeam_genambig)[1]
loci_all  <- Loci(hornbeam_genambig)
samps_p1  <- Samples(hornbeam_genambig, populations = pop1)

cat("First population:", pop1, "\n")
cat("All loci:", paste(loci_all, collapse = ", "), "\n\n")

# Count distinct non-null alleles per locus for pop1 (excluding missing)
allele_counts <- sapply(loci_all, function(L) {
  s_ok <- samps_p1[!isMissing(hornbeam_genambig, samps_p1, L)]
  if (length(s_ok) == 0L) return(Inf)
  all_a <- unlist(lapply(s_ok, function(s) Genotype(hornbeam_genambig, s, L)))
  length(unique(all_a[all_a != 0]))  # 0 is sometimes used for missing
})
cat("Distinct alleles per locus in", pop1, ":\n")
print(allele_counts)

best_locus <- names(which.min(allele_counts))
na_best    <- allele_counts[best_locus]
cat("\nChosen locus:", best_locus, "(", na_best, "distinct alleles )\n\n")

# --- Subset object to pop1 + best_locus ------------------------------------
sub_obj <- hornbeam_genambig[samps_p1, best_locus]

# Show genotype space size
m2  <- unique(as.vector(Ploidies(sub_obj)))
na1 <- na_best + 1L
ng  <- choose(na1 + m2 - 1L, m2)
cat(sprintf("Ploidy m2 = %d,  na = %d,  ng = %.0f\n", m2, na_best, ng))
cat(sprintf("Code path: %s\n\n", if (ng > 50000) "SPARSE" else "DENSE"))

# --- Run with profiling (time + peak memory) --------------------------------
# Memory before
gc(reset = TRUE)
mem_before <- sum(gc()[, 2])   # sum of used Vcells (MB) from gc

t0 <- proc.time()

result <- de_silva_freq_EM(sub_obj, self = 0.5, loci = best_locus)

elapsed <- proc.time() - t0

# Memory after
gc_after   <- gc()
mem_after  <- sum(gc_after[, 6])   # max used (Vcells) in Mb column 6

cat("=== Timing ===\n")
cat(sprintf("  elapsed : %.2f s\n", elapsed["elapsed"]))
cat(sprintf("  user    : %.2f s\n", elapsed["user.self"]))
cat(sprintf("  sys     : %.2f s\n", elapsed["sys.self"]))

cat("\n=== Memory (from gc()) ===\n")
print(gc_after)

cat("\n=== Result ===\n")
freq_cols <- grep(paste0("^", best_locus, "\\."), names(result))
cat(sprintf("Allele frequency columns: %d\n", length(freq_cols)))
cat(sprintf("Freq sum for %s: %.8f\n", pop1, sum(result[pop1, freq_cols])))
print(round(result[pop1, freq_cols], 4))
