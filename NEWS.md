# dng NEWS

## dng 1.0.0

- Added `gsplitn()` for split-normal CDF and log-density gradients.
- Added and expanded `testthat` coverage for split-normal and split-t
  distributions, moments, random generation, and gradients.
- Fixed split-t density and kurtosis bugs.
- Modernized package metadata and citation format.
- Kept R documentation generated from C++ roxygen comments.

## dng 0.2.1

- Added split-t gradient support through `gsplitt()`.
- Improved C++ argument recycling with `rep_len()`.
- Fixed `qsplitn()` issues.
- Consolidated distribution and moment documentation.
- Refreshed package metadata, generated Rcpp exports, and tests.

## dng 0.2.0

- Added split-normal distribution functions: `dsplitn()`, `psplitn()`,
  `qsplitn()`, and `rsplitn()`.
- Added split-normal moment functions for mean, variance, skewness, and
  kurtosis.
- Added split-normal documentation and examples.

## dng 0.1.1

- Added split-t distribution functions: `dsplitt()`, `psplitt()`, `qsplitt()`,
  and `rsplitt()`.
- Added split-t moment functions for mean, variance, skewness, and kurtosis.
- Added initial package documentation, citation metadata, and package
  infrastructure.
- Improved Rcpp compatibility and performance.
