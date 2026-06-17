## Submission

This is a new release (first submission) of rolescry 0.1.0.

## Test environments

* Local: Windows 11, R 4.5.3, `R CMD check --as-cran` (with manual and vignettes).
* win-builder: R-release (4.6.0) and R-devel (r90166, 2026-06-16).
* R-hub v2 (GitHub Actions): linux (R-devel), windows (R-devel), clang-asan, valgrind — all OK.
* R-universe: Windows, macOS, and Linux builds — OK.

## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE is from the CRAN incoming feasibility check:

* New submission.
* Possibly misspelled words in DESCRIPTION: "Boynukara", "biostatistics",
  "inspice", "nomen". These are all correct: "Boynukara" is the author's
  surname; "biostatistics" is the standard domain term; "inspice" and "nomen"
  are Latin words from the project motto ("Data inspice, non nomen").

There are no other notes on r-release or r-devel. The PDF and HTML manuals,
vignette, examples, and tests all build and pass.

## Reverse dependencies

This is a new package, so there are no reverse dependencies.
