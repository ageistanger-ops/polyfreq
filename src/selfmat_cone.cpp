// selfmat_cone.cpp
//
// Rcpp implementation of the sparse selfing matrix for the union-phenotype cone.
// Replaces .selfmat_cone_batch() in R/helpers_desilva.R.
//
// All pre-computation (gamete positions, unordered pair indices, hash weights,
// cone hashes) is performed in R and passed in as arguments.  The C++ function
// handles only the heavy inner loops: gamete extraction, inline merge-sort,
// hash lookup, and triplet accumulation.
//
// OpenMP parallelises the chunk loop when nthreads > 1.  Each thread owns a
// local triplet vector; results are merged after the parallel section.
//
// Reference: De Silva HN et al. (2005) Heredity 95:327-334.

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Helper: compare two Triplets for sorting by (par, off)
// ---------------------------------------------------------------------------
struct Triplet {
  int    par;
  int    off;
  double wt;
};

static inline bool triplet_less(const Triplet& a, const Triplet& b) {
  if (a.par != b.par) return a.par < b.par;
  return a.off < b.off;
}

// ---------------------------------------------------------------------------
// selfmat_cone_cpp
//
// Arguments (all pre-computed in R, 0-based indices where noted):
//   ag_cone       IntegerMatrix (n_cone x m2)  — cone genotypes (1-based allele pos)
//   gam_pos       IntegerMatrix (m x n_gam)    — gamete position combos (1-based col idx)
//   idx1          IntegerVector (n_pairs)       — first gamete index, 0-based
//   idx2          IntegerVector (n_pairs)       — second gamete index, 0-based
//   pair_wt       IntegerVector (n_pairs)       — 1 (homo) or 2 (hetero)
//   hash_weights  NumericVector (m2)            — polynomial hash weights
//   cone_hashes   NumericVector (n_cone)        — hash of each cone row
//   m             int                           — half-ploidy (m2/2)
//   chunk_size    int                           — parents per chunk
//   nthreads      int                           — OpenMP threads (1 = serial)
//
// Returns a List with integer vectors i, j and numeric vector x (1-based)
// suitable for Matrix::sparseMatrix(i, j, x, dims = c(n_cone, n_cone)).
// ---------------------------------------------------------------------------

// [[Rcpp::export(.selfmat_cone_cpp)]]
List selfmat_cone_cpp(IntegerMatrix ag_cone,
                      IntegerMatrix gam_pos,
                      IntegerVector idx1,
                      IntegerVector idx2,
                      IntegerVector pair_wt,
                      NumericVector hash_weights,
                      NumericVector cone_hashes,
                      int m,
                      int chunk_size,
                      int nthreads) {

  const int n_cone  = ag_cone.nrow();
  const int m2      = 2 * m;
  const int n_gam   = gam_pos.ncol();
  const int n_pairs = idx1.size();
  const int n_chunks = (n_cone + chunk_size - 1) / chunk_size;

  // --- Build cone hash table (hash -> 1-based cone row index) ---------------
  std::unordered_map<double, int> hash_table;
  hash_table.reserve(static_cast<size_t>(n_cone) * 2);
  for (int g = 0; g < n_cone; g++) {
    hash_table[cone_hashes[g]] = g + 1;   // 1-based
  }

  // --- Thread-local triplet storage -----------------------------------------
  // Each OpenMP thread accumulates its own vector; merged after parallel loop.
  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (nthreads > 0) ? nthreads : omp_get_max_threads();
#endif
  std::vector< std::vector<Triplet> > thread_triplets(actual_threads);
  for (int t = 0; t < actual_threads; t++) {
    // Reserve a rough upper bound to avoid repeated reallocations.
    thread_triplets[t].reserve(
      static_cast<size_t>(chunk_size) * n_pairs / 4
    );
  }

  // --- Parallel chunk loop --------------------------------------------------
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 1) num_threads(actual_threads)
#endif
  for (int ck = 0; ck < n_chunks; ck++) {

    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif

    const int cs = ck * chunk_size;
    const int ce = std::min(cs + chunk_size, n_cone);
    const int B  = ce - cs;

    // --- Step 1: build gamete matrix gmat[B * n_gam][m] (row-major) ---------
    // gmat[(b * n_gam + ig) * m + k] = allele at position gam_pos(k, ig)-1
    //                                   of parent ag_cone[cs + b, ...]
    std::vector<int> gmat(static_cast<size_t>(B) * n_gam * m);

    for (int b = 0; b < B; b++) {
      const int parent_row = cs + b;
      for (int ig = 0; ig < n_gam; ig++) {
        const int gmat_base = (b * n_gam + ig) * m;
        for (int k = 0; k < m; k++) {
          // gam_pos is 1-based column index into ag_cone
          const int col = gam_pos(k, ig) - 1;   // 0-based column of ag_cone
          gmat[gmat_base + k] = ag_cone(parent_row, col);
        }
      }
    }

    // --- Steps 2-4: for each parent b and each unordered pair (ig1, ig2) ---
    //   - inline merge-sort the two m-element gametes into off[m2]
    //   - hash off[] and look up in hash_table
    //   - if found, push Triplet to thread-local vector

    std::vector<Triplet>& local_trips = thread_triplets[tid];

    for (int b = 0; b < B; b++) {
      const int global_par = cs + b + 1;   // 1-based
      const int gam_base   = b * n_gam;

      for (int p = 0; p < n_pairs; p++) {
        const int ig1      = idx1[p];   // 0-based gamete index
        const int ig2      = idx2[p];
        const int* g1      = &gmat[(gam_base + ig1) * m];
        const int* g2      = &gmat[(gam_base + ig2) * m];

        // Inline merge of two sorted m-element arrays into off[m2]
        int off[8];   // m2 <= 8 for octoploid; stack-allocated
        {
          int ia = 0, ib = 0;
          for (int k = 0; k < m2; k++) {
            if (ia < m && (ib >= m || g1[ia] <= g2[ib])) {
              off[k] = g1[ia++];
            } else {
              off[k] = g2[ib++];
            }
          }
        }

        // Polynomial hash
        double h = 0.0;
        for (int k = 0; k < m2; k++) {
          h += off[k] * hash_weights[k];
        }

        // Lookup
        auto it = hash_table.find(h);
        if (it != hash_table.end()) {
          local_trips.push_back({global_par, it->second,
                                 static_cast<double>(pair_wt[p])});
        }
      }
    }

    // --- Aggregate triplets within this chunk for this thread ---------------
    // Sort the newly added triplets (for this chunk) by (par, off), then
    // aggregate runs with the same key by summing weights.
    //
    // We track the start position before the chunk was added:
    // find the boundary between previously accumulated and new triplets.
    // Actually, we sort and aggregate PER CHUNK to keep memory bounded.
    // We use a simple approach: mark the start, sort the suffix, aggregate.
    //
    // To do this correctly without marking: accumulate raw into a temporary
    // vector per chunk, then merge into local_trips after aggregation.
    // We re-do this by using the raw push_back approach and doing a final
    // sort+aggregate after the parallel section.
    // (Sorting per-chunk would require knowing the start offset per chunk,
    // which is complex with OpenMP.  The final sort is O(N log N) total,
    // acceptable since total triplets << B * n_pairs.)
  }

  // --- Merge thread-local triplets and aggregate ----------------------------
  // Compute total size, then merge all thread vectors into one.
  size_t total = 0;
  for (int t = 0; t < actual_threads; t++)
    total += thread_triplets[t].size();

  std::vector<Triplet> all_trips;
  all_trips.reserve(total);
  for (int t = 0; t < actual_threads; t++) {
    all_trips.insert(all_trips.end(),
                     thread_triplets[t].begin(),
                     thread_triplets[t].end());
    // Free thread-local memory immediately
    std::vector<Triplet>().swap(thread_triplets[t]);
  }

  // Sort by (par, off)
  std::sort(all_trips.begin(), all_trips.end(), triplet_less);

  // Aggregate runs with same (par, off) into std::vectors, then wrap for R
  std::vector<int>    si_out, sj_out;
  std::vector<double> sx_out;
  si_out.reserve(all_trips.size());
  sj_out.reserve(all_trips.size());
  sx_out.reserve(all_trips.size());

  for (size_t t = 0; t < all_trips.size(); ) {
    size_t t2    = t + 1;
    double sum_w = all_trips[t].wt;
    while (t2 < all_trips.size() &&
           all_trips[t2].par == all_trips[t].par &&
           all_trips[t2].off == all_trips[t].off) {
      sum_w += all_trips[t2].wt;
      t2++;
    }
    si_out.push_back(all_trips[t].par);
    sj_out.push_back(all_trips[t].off);
    sx_out.push_back(sum_w);
    t = t2;
  }

  return List::create(Named("i") = IntegerVector(si_out.begin(), si_out.end()),
                      Named("j") = IntegerVector(sj_out.begin(), sj_out.end()),
                      Named("x") = NumericVector(sx_out.begin(), sx_out.end()));
}
