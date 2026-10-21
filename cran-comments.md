## Release summary

rolescry 0.3.0 (update of 0.2.0, on CRAN since 2026-07-19).

Changes: three derived decision rules (group variable by normalized entropy,
paired vs agreement measurements by an MDL shift code, ordinal scale items by
an MDL one-factor code); two new exported functions, `verify_roles()` (tests a
declared role assignment against the data) and `harmonize_roles()` (aligns
detected roles across several data frames). See NEWS.md.

## Test environments

* local: Windows 11, R 4.6.1 -- 0 errors, 0 warnings, 0 notes
* win-builder: R-devel and R-release -- PASTE_STATUS
* macOS builder (R-release) -- PASTE_STATUS

## R CMD check results

0 errors | 0 warnings | 0 notes locally.

The only NOTE on the remote builders is the CRAN incoming feasibility check
listing possibly-misspelled words in DESCRIPTION: "Boynukara" (author's
surname), "biostatistics" (domain term), and "inspice"/"nomen" (Latin, from the
motto "Data inspice, non nomen"). All are correct.

## Reverse dependencies

There are no reverse dependencies on CRAN.
