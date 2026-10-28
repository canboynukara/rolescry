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

* local: Windows 11, R 4.5.3 -- 0 errors, 0 warnings, 1 note ("unable to
  verify current time", environment only); Linux x86_64, R 4.3.3 -- tests OK
* win-builder: R-release 4.6.1 -- Status: OK
* win-builder: R-devel -- submitted 2026-10-28 (result pending at the time of writing)
* macOS builder: r-release-macosx-arm64, R 4.6.1, Apple M1, clang-1700 --
  Status: OK (errors:no, warnings:no, notes:no) -- the platform of the
  reported failure

## R CMD check results

0 errors | 0 warnings | 0 notes (win-builder, macOS arm64)

## Reverse dependencies

There are no reverse dependencies on CRAN.
