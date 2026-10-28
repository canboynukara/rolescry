## Release summary

rolescry 0.3.1 -- a fix release requested by the CRAN team (email of
2026-10-28, deadline 2026-11-13): the tests of 0.3.0 failed on
r-release-macos-arm64 (M1mac) and in the noLD additional check.

Cause: in one test fixture two candidate columns had exactly symmetric
evidence (the same 2x2 statistic computed two ways), and the canonical pick
compared the floating-point scores exactly, so a last-ulp difference on
arm64 / without long doubles decided the pick instead of the deterministic
value-fingerprint tie-break.

Fix: `.canonical_pick()` now treats scores within 64 machine epsilons as a
tie and resolves it by the value fingerprint on every platform; the fixture
was changed so that it is genuinely decisive. No change to any detection
rule. See NEWS.md.

## Test environments

* local: Linux (x86_64), R 4.3.3 -- 0 errors, 0 warnings, 0 notes
* win-builder: R-release and R-devel -- (fill in after submission)
* macOS builder (R-release, arm64) -- (fill in after submission)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
