# helper.R — shared test fixtures and utilities
# Loaded automatically by testthat before all test files.

library(polysat)

# Load the simgen test dataset (300 tetraploid samples, 3 loci, 3 populations)
simgen_path <- system.file("testdata", "simgen.Rdata", package = "polyfreq")
if (nzchar(simgen_path)) {
  load(simgen_path)   # loads 'simgen' into the test environment
} else {
  stop("Cannot find simgen.Rdata in inst/testdata/ — is the package installed?")
}

# Compute the polysat gold-standard result once for the whole test session
# to avoid repeating the (slow) C-based computation in every test.
POLYSAT_RESULT <- polysat::deSilvaFreq(simgen, self = 0.1)
