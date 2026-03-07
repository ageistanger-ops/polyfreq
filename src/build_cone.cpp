// build_cone.cpp
//
// Rcpp implementation of the union-phenotype cone construction.
// Replaces the three heavy pure-R loops inside .build_cone() in
// R/de_silva_freq_EM.R:
//
//   1. Per-phenotype GENLIST_r() calls (one per unique observed phenotype)
//   2. do.call(rbind, ...) + !duplicated() global deduplication
//   3. n_cone-iteration loop computing indexf() for every cone genotype
//   4. sort/unique/match for compact phenotype mapping
//
// The polysat S4 accessor calls (Genotype, Samples, ...) that collect
// unique_phenos remain in R because they require the polysat class system.
//
// Reference: De Silva HN et al. (2005) Heredity 95:327-334.

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// choose_int: exact C(n, k) for small non-negative integers.
// Returns 0 for any out-of-range input (n < 0, k < 0, or k > n).
// ---------------------------------------------------------------------------
static inline int choose_int(int n, int k) {
  if (n < 0 || k < 0 || k > n) return 0;
  if (k == 0 || k == n)        return 1;
  if (k > n - k) k = n - k;     // use min(k, n-k) to keep intermediate small
  long long r = 1;
  for (int i = 0; i < k; ++i)
    r = r * (n - i) / (i + 1);
  return static_cast<int>(r);
}

// G_r(q, n) = C(n + q, q + 1)  — de Silva G function
static inline int G_r(int q, int n) {
  return choose_int(n + q, q + 1);
}

// ---------------------------------------------------------------------------
// indexf_cpp
//
// Maps a sorted non-null allele phenotype to its INDEXF integer (1-based).
//
// Arguments:
//   af1   pointer to sorted 1-based allele positions (values 1..na)
//   m     number of non-null alleles (0 for null phenotype)
//   na    total number of non-null alleles at this locus
//
// Implements .make_indexf(m2) from de_silva_freq_EM.R, which is a direct
// R translation of the INDEXF subroutine from de Silva et al. (2005).
// All af1[i] references are 0-based (C++) equivalents of the 1-based R code.
// ---------------------------------------------------------------------------
static inline int indexf_cpp(const int* af1, int m, int na) {
  int x = 1;
  if (m == 0) return x;
  if (m == 1) {
    x += af1[0];          // R: x <- 1 + af1[1L]  (1-based)
    return x;
  }
  // m >= 2
  for (int q = 1; q <= m - 1; ++q)
    x += G_r(q - 1, na - q + 1);

  x += G_r(m - 1, na + 1 - m) - G_r(m - 1, na + 2 - m - af1[0]);

  if (m > 2) {
    // R: for q in 1:(m-2): G_r(q, na-q-af1[m-q-1]) - G_r(q, na+1-q-af1[m-q])
    // 1-based R indices af1[m-q-1] and af1[m-q]
    // → 0-based C++:   af1[m-q-2]             and af1[m-q-1]
    for (int q = 1; q <= m - 2; ++q)
      x += G_r(q, na - q - af1[m - q - 2]) -
           G_r(q, na + 1 - q - af1[m - q - 1]);
  }

  // R: x <- x + af1[m] - af1[m-1L]  (1-based)
  // → 0-based C++: af1[m-1] - af1[m-2]
  x += af1[m - 1] - af1[m - 2];
  return x;
}

// ---------------------------------------------------------------------------
// build_cone_cpp
//
// Arguments:
//   unique_phenos  R List of IntegerVector — one vector per unique observed
//                  phenotype; each vector holds sorted 1-based positions in
//                  alleles[] (values 1..na, no null, length k >= 0).
//   na             int — number of non-null alleles at this locus
//   na1            int — na + 1 (null sentinel index in local 1..na1 space)
//   m2             int — ploidy level (must be even)
//
// Returns a List:
//   ag_cone        IntegerMatrix (n_cone x m2)  — cone genotypes; each row
//                  is a sorted non-decreasing sequence of global 1-based
//                  allele-position indices (values in 1..na1, na1 = null)
//   compact_ph_idx IntegerVector (n_cone)        — compact phenotype index
//                  (1..n_uph) for each cone genotype
//   unique_ph      IntegerVector (n_uph)         — raw INDEXF value for each
//                  compact phenotype group
// ---------------------------------------------------------------------------
// [[Rcpp::export(.build_cone_cpp)]]
List build_cone_cpp(List   unique_phenos,
                    int    na,
                    int    na1,
                    int    m2) {

  const int n_phenos = unique_phenos.size();

  // Polynomial hash for deduplication.
  // base = na1 + 1 > na1 (max allele position) so that
  //   h = sum_j(ag[j] * base^j) is injective over valid sorted sequences.
  const double hash_base = static_cast<double>(na1 + 1);
  std::vector<double> hw(m2);
  hw[0] = 1.0;
  for (int j = 1; j < m2; ++j) hw[j] = hw[j - 1] * hash_base;

  // Flat accumulators for all genotype rows (before deduplication).
  // all_rows: stride-m2, global 1-based allele positions.
  std::vector<int>    all_rows;
  std::vector<double> all_hashes;
  all_rows.reserve(60000 * m2);
  all_hashes.reserve(60000);

  std::vector<int> ag1(m2);    // scratch: current local genotype (1-based)

  // -------------------------------------------------------------------------
  // Step 1: for each observed phenotype, enumerate all m2-multisets from
  //         { pheno[0], ..., pheno[k-1], null } using inline GENLIST,
  //         map local indices to global, and push to accumulators.
  // -------------------------------------------------------------------------
  for (int pi = 0; pi < n_phenos; ++pi) {

    const IntegerVector pheno = unique_phenos[pi];  // 1-based positions 1..na
    const int k               = pheno.size();
    const int local_na1       = k + 1;

    // Local → global mapping (0-based indexing into gmap):
    //   local index v (1-based, 1..local_na1) → gmap[v - 1]
    //   v = 1..k   → pheno[v-1]  (1-based global position in 1..na)
    //   v = k+1    → na1         (null sentinel)
    // gmap is monotonically non-decreasing (pheno is sorted, na1 > all pheno
    // values), so global genotype rows inherit the non-decreasing order of
    // GENLIST — required for the INDEXF computation below.
    std::vector<int> gmap(local_na1);
    for (int j = 0; j < k; ++j) gmap[j] = pheno[j];
    gmap[k] = na1;

    // Inline GENLIST: enumerate all m2-multisets from {1..local_na1}.
    // Iterative algorithm identical to GENLIST_r in R/helpers_desilva.R.
    for (int j = 0; j < m2; ++j) ag1[j] = 1;

    // Emit helper: map ag1 to global, compute hash, push to accumulators.
    auto emit = [&]() {
      double h = 0.0;
      for (int j = 0; j < m2; ++j) {
        const int gv = gmap[ag1[j] - 1];
        all_rows.push_back(gv);
        h += static_cast<double>(gv) * hw[j];
      }
      all_hashes.push_back(h);
    };

    emit();                   // first genotype: (1, 1, ..., 1) in local space
    int a = m2 - 1;           // rightmost position (0-based)
    while (a >= 0) {
      if (ag1[a] == local_na1) {
        a--;
      } else {
        ag1[a]++;
        for (int j = a + 1; j < m2; ++j) ag1[j] = ag1[a];
        emit();
        a = m2 - 1;
      }
    }
  }

  const int n_total = static_cast<int>(all_hashes.size());

  // -------------------------------------------------------------------------
  // Step 2: deduplicate by hash, keeping the first occurrence of each unique
  //         genotype across all phenotype closures.
  // -------------------------------------------------------------------------
  std::unordered_map<double, int> seen;
  seen.reserve(static_cast<size_t>(n_total) * 2);

  std::vector<int> keep_idx;
  keep_idx.reserve(n_total / 2);

  for (int i = 0; i < n_total; ++i) {
    if (seen.emplace(all_hashes[i], static_cast<int>(keep_idx.size())).second)
      keep_idx.push_back(i);
  }

  const int n_cone = static_cast<int>(keep_idx.size());

  // Build ag_cone matrix from kept rows.
  IntegerMatrix ag_cone(n_cone, m2);
  for (int g = 0; g < n_cone; ++g) {
    const int base = keep_idx[g] * m2;
    for (int j = 0; j < m2; ++j)
      ag_cone(g, j) = all_rows[base + j];
  }

  // Free accumulators now that ag_cone is built.
  { std::vector<int>().swap(all_rows); }
  { std::vector<double>().swap(all_hashes); }

  // -------------------------------------------------------------------------
  // Step 3: compute INDEXF for each cone genotype.
  //
  // Each ag_cone row is sorted non-decreasing (guaranteed by GENLIST +
  // monotone gmap).  Extract unique non-null alleles (values < na1) from
  // each row — these are already sorted and deduplicated by construction.
  // -------------------------------------------------------------------------
  std::vector<int> raw_ph(n_cone);
  std::vector<int> af1_buf(m2);

  for (int g = 0; g < n_cone; ++g) {
    int naf1 = 0;
    if (ag_cone(g, 0) != na1) {
      af1_buf[0] = ag_cone(g, 0);
      naf1 = 1;
      for (int j = 1; j < m2; ++j) {
        const int v = ag_cone(g, j);
        if (v == na1) break;
        if (v > af1_buf[naf1 - 1]) af1_buf[naf1++] = v;
      }
    }
    raw_ph[g] = indexf_cpp(af1_buf.data(), naf1, na);
  }

  // -------------------------------------------------------------------------
  // Step 4: compact phenotype mapping.
  // Sort unique raw INDEXF values → 1-based compact indices 1..n_uph.
  // -------------------------------------------------------------------------
  std::vector<int> uph = raw_ph;
  std::sort(uph.begin(), uph.end());
  uph.erase(std::unique(uph.begin(), uph.end()), uph.end());
  const int n_uph = static_cast<int>(uph.size());

  std::unordered_map<int, int> ph_map;
  ph_map.reserve(static_cast<size_t>(n_uph) * 2);
  for (int i = 0; i < n_uph; ++i) ph_map[uph[i]] = i + 1;   // 1-based

  IntegerVector compact_ph_idx(n_cone);
  for (int g = 0; g < n_cone; ++g)
    compact_ph_idx[g] = ph_map[raw_ph[g]];

  IntegerVector unique_ph(uph.begin(), uph.end());

  return List::create(
    Named("ag_cone")        = ag_cone,
    Named("compact_ph_idx") = compact_ph_idx,
    Named("unique_ph")      = unique_ph
  );
}
