## Update

This is an update of rolescry, 0.1.0 -> 0.2.0.

It ships (1) a correctness fix to the name-blindness contract (column typing no
longer consults column names when `name_bonus = NULL`), (2) more robust survival
event-code normalization, and (3) an information-theoretic re-founding of the
detection engine. The name-blindness fix is a behaviour change and is documented
in NEWS.md. Version 0.1.0 reached CRAN roughly one month ago; this update is
motivated by the correctness fix.

## Test environments

* Local: Windows 11, R 4.5.3, `R CMD check --as-cran`  -- MAINTAINER TO CONFIRM.
* win-builder: R-release and R-devel  -- MAINTAINER TO RUN AND CONFIRM before submission.
* R-hub v2 (GitHub Actions): linux (R-devel), windows (R-devel)  -- MAINTAINER TO CONFIRM.

## R CMD check results

Expected: 0 errors | 0 warnings | 1 note.

The single NOTE is the CRAN incoming feasibility check listing possibly-misspelled
words in DESCRIPTION: "Boynukara" (the author's surname), "biostatistics" (the
standard domain term), and "inspice"/"nomen" (Latin, from the motto
"Data inspice, non nomen"). All are correct.

## Reverse dependencies

There are no reverse dependencies on CRAN.
