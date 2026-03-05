# polyfreq — Development Notes

## Project Goal
Pure-R reimplementation of the de Silva et al. (2005) EM algorithm for allele frequency estimation in polyploids, originally implemented with C++/Rcpp helpers in the `polysat` package.

## Success Criteria
1. De Silva frequency estimation fully transferred to pure R (no C/C++ code)
2. R vectorization and fast idioms applied wherever possible
3. Numerically identical results to `polysat::deSilvaFreq` (tolerance: 1e-8)
4. Runtime equal to or faster than the polysat implementation

## Architecture
- Uses `polysat::genambig` class and accessors — do NOT reimplement
- Dependencies: polysat (class system + simpleFreq), cli (errors), Matrix (kept for future use)
- 5 exported helper functions in `R/helpers_desilva.R`: G_r, INDEXG_r, GENLIST_r, RANMUL_r, SELFMAT_r
- 1 exported main function + 3 internal helpers in `R/de_silva_freq_EM.R`
- 106 pre-written tests in `tests/testthat/`

## Key Performance Strategies
- `G_r`: Use `choose()` (vectorized C call) instead of loop
- `RANMUL_r`: Vectorize with `tabulate()` and `lfactorial()`
- `SELFMAT_r` (bottleneck): Batch INDEXG via vectorized `choose()`, `combn()` pre-computation, `tabulate()` for row fills
- EM loop: Vectorized `rvec` via `pa[ag]` matrix indexing; start with per-iteration `solve()`, optimize to pre-computed inverse if needed

## R Installation
`C:\Program Files\R\R-4.4.3\bin`

## Reference Materials
- `polyfreq_material/6800728.pdf` — Original de Silva et al. (2005) publication
- `polyfreq_material/polysat.pdf` — polysat package manual
- `polyfreq_material/R/` — polysat R source code (esp. population_stats.R lines 110-422)
- `polyfreq_material/src/deSilvaFuncs.cpp` — C++ helpers to reimplement
