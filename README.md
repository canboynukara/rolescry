# rolescry

> Name-blind variable-role detection by data signature. *Data inspice, non nomen* -- inspect the data, not the name.

`rolescry` assigns statistical **roles** to the columns of a tabular dataset --
group variable, continuous/binary outcome, survival time and event, paired and
agreement measurement pairs, repeated measures, scale items, subject identifier,
and covariates -- using only each column's **information-theoretic signature**
(Shannon entropy, normalized mutual information, distributional shape and
inter-column structure), never the column names. Renaming every column to
`col_1, col_2, ...` does not change the result. No large language models, no
external data transmission; detection is deterministic.

Detection is backed by two structural invariance guarantees -- renaming
(RELABEL) and reordering (S_n) columns never change a result -- and a single
pre-registered, held-out confirmatory run (OSF osf.io/8ecau) validated the
estimator on a synthetic data-generating process. Real-data and external
validity are a separate, ongoing question.

## Relationship to MDStatR

`rolescry` began as a variable-role detection component written during development of
MDStatR, an unreleased personal research codebase by the same author. It was separated
into an independent project, and for version 0.2.0 its detection engine was rebuilt on
an information-theoretic basis: minimum description length scoring over a
distance-correlation dependency backbone.

For anyone evaluating or using `rolescry`, the practical position is:

- **No dependency.** `rolescry` requires only base R (`stats`, `utils`). It does not
  import, link to or call MDStatR, and MDStatR is not needed to install, run, test or
  understand it.
- **No overlapping distributed software.** MDStatR is not publicly released and is not
  installable; it is archived on Zenodo as a record only. No other package duplicates
  `rolescry`'s functionality.
- **Its own contracts and evidence.** The relabelling-invariance and permutation-
  equivariance guarantees, the test suite that enforces them, and the pre-registered
  evaluation belong to `rolescry` and have no counterpart in MDStatR.

The MDStatR reference is retained to attribute intellectual origin, nothing more.

## Installation

From CRAN:

```r
install.packages("rolescry")
```

From [r-universe](https://r-universe.dev) (development builds):

```r
install.packages("rolescry", repos = "https://canboynukara.r-universe.dev")
```

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("canboynukara/rolescry")
```

The package needs only **base R + `stats`**. Optional packages
(`readxl`/`openxlsx`/`haven` for file reading; `moments`/`diptest`/`stringdist`
for extra refinements) are used only if installed.

## Quick start

```r
library(rolescry)

set.seed(1)
d <- data.frame(
  arm  = rep(c(0, 1), each = 50),   # group
  pre  = rnorm(100, 10, 2),         # paired with post
  post = rnorm(100, 11, 2),
  resp = rbinom(100, 1, 0.4)        # binary outcome
)

res <- detect_roles(d)
res
res$roles$group_var$columns
summary(res)
```

### The name-blindness guarantee

Detection is purely mathematical by default (`name_bonus = NULL`):

```r
pos <- function(res, dat) match(res$roles$paired_pairs$columns, names(dat))
d_blind <- setNames(d, paste0("col_", seq_along(d)))
identical(pos(detect_roles(d), d), pos(detect_roles(d_blind), d_blind))
#> TRUE  -- the SAME columns (by position) are detected, named or col_N
```

Column names can be used only as a small, **capped** tie-breaker (at most a
+10 point nudge, i.e. <= 10%) by passing a keyword dictionary; the mathematical
signature still dominates:

```r
detect_roles(d, name_bonus = rolescry_default_name_bonus())
```

### Header-aware loading

```r
df <- read_data("messy_export.xlsx")   # auto-detects the header row
```

## How it works

`detect_roles()` types each column from its values (`.build_var_info`), scores
candidate roles with information-theoretic and distributional signatures
(`compute_nmi()` exposes the normalized mutual information directly), and returns
a structured `role_detection` object with per-role confidence and a component
breakdown. See `vignette("rolescry")` for the method and the name-blind
guarantee.

## Citation & attribution

`rolescry` has its own archival DOI: **10.5281/zenodo.21003941** (Zenodo). It is
derived from the MDStatR engine: Boynukara, C. (2026). *MDStatR (v2.1.0 Veritas).*
Zenodo. https://doi.org/10.5281/zenodo.20707791

Run `citation("rolescry")` to cite the package (with its archival DOI) and its
parent engine.

## License

Apache License 2.0.
