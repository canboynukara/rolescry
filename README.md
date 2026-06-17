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

Extracted from the MDStatR biostatistics engine.

## Installation

From [r-universe](https://r-universe.dev):

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

Derived from Boynukara, C. (2026). MDStatR (v2.1.0 Veritas). Zenodo. https://doi.org/10.5281/zenodo.20707791

Run `citation("rolescry")` to cite the package and its parent engine.

## License

Apache License 2.0 (inherited from the parent MDStatR project). See `LICENSE.md`.
