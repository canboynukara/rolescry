# rolescry

> Name-blind variable-role detection by data signature. *Data inspice, non nomen* -- inspect the data, not the name.

`rolescry` assigns statistical **roles** to the columns of a tabular dataset
from **any domain** -- a clinical trial, a loan book, a firm-month panel, a
survey --
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
into an independent project and is released on its own terms.

Version 0.2.0 rebuilt part of the detection engine on an information-theoretic basis
(measurement-pair and outcome detection on minimum description length scoring over a
distance-correlation backbone, an additive-noise-model direction test, a bit-margin
abstention rule). Version 0.3.0 moves three more roles onto derived criteria: `group_var`
(normalized Shannon entropy), paired versus agreement pairs (the sign of the MDL shift bits)
and ordinal `scale_items` (a one-factor Gaussian two-part code). Three rules remain on the
original ported signature scorers -- survival, wide `repeated_measures` and binary-item
`scale_items` -- because derived replacements were measured and rejected under a
no-regression rule (see Known limitations); the source files record provenance individually.

For anyone evaluating or using `rolescry`, the practical position is:

- **No dependency.** `rolescry` requires only base R (`stats`, `utils`). It does not
  import, link to or call MDStatR, and MDStatR is not needed to install, run, test or
  understand it.
- **No overlapping distributed software.** MDStatR is not publicly released and is not
  installable; it is archived on Zenodo as a record only. No other package duplicates
  `rolescry`'s functionality.
- **Its own contracts and evidence.** The relabeling-invariance and permutation-
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
n <- 200
pre <- rnorm(n, 10, 2)
d <- data.frame(
  arm  = rep(c("control", "treatment"), each = n / 2),  # balanced 2-level grouping
  pre  = pre,                                           # a measurement ...
  post = pre + rnorm(n, 1, 1),                          # ... repeated after (paired)
  z1   = rnorm(n),                                      # two standardized predictors
  z2   = rnorm(n)
)
d$resp <- rbinom(n, 1, plogis(1.2 * d$z1 + 1.2 * d$z2)) # binary outcome caused by z1, z2

res <- detect_roles(d)
res$roles$group_var$columns        #> "arm"
res$roles$paired_pairs$columns     #> "pre" "post"
res$roles$outcome_binary$columns   #> "resp"
res$roles$covariate$columns        #> "z1" "z2"   (columns no other role claimed)
summary(res)
```

### The name-blindness guarantee

Detection is purely mathematical by default (`name_bonus = NULL`):

```r
pos <- function(res, dat) match(res$roles$outcome_binary$columns, names(dat))
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

### Beyond biostatistics

The detectors are information-theoretic, not domain-bound. The same calls work on
financial data:

```r
# a loan book: duration in months, rare default indicator
set.seed(31); n <- 1000
dur <- round(rweibull(n, 1.2, 36)) + 1
def <- rbinom(n, 1, plogis(-1 - 0.08 * dur))
loans <- data.frame(amount = rlnorm(n, 10, 0.6), rate = rnorm(n, 0.08, 0.02),
                    dur = dur, def = def)
detect_roles(loans)$roles$time_variable$columns   #> "dur"   (duration/event pair)

# a daily price series: the Date column is a calendar axis (time_index)
prices <- data.frame(dt = as.Date("2024-01-01") + 0:199,
                     price = cumsum(rnorm(200)) + 100)
detect_roles(prices)$roles$time_index$columns     #> "dt"

# finance keyword dictionary for the optional name hint
detect_roles(loans, name_bonus = rolescry_default_name_bonus("finance"))
```

## Verifying a declared analysis plan

`verify_roles(data, claim)` turns the detector around: you state which columns
play which role, and the data are asked whether they can support it. The check is
name-blind (the claim is mapped to column positions once) and returns one verdict
per claimed role.

```r
set.seed(21); n <- 300
arm <- rep(c("drug", "placebo"), each = n / 2)
base <- rnorm(n, 150, 15)
rct <- data.frame(id = sample.int(99999, n), arm, age = round(rnorm(n, 60, 10)),
                  sbp_base = base, sbp_end = base - 10 * (arm == "drug") + rnorm(n, 0, 8))
rct$y <- rbinom(n, 1, plogis(-1 + 0.05 * (rct$sbp_end - 140)))
verify_roles(rct, list(subject_id = "id", group_var = "arm", paired_pairs = c("age", "sbp_end"),
                       outcome_binary = "y"))
#> <role_verification> 4 claims: 2 testable, 1 contradicted, 1 consistent, 2 unverifiable
#>            role      columns      verdict      sub_reason   statistic
#>      subject_id           id   consistent            <NA> unique=TRUE
#>       group_var          arm unverifiable identifiability     B=1.000
#>    paired_pairs age, sbp_end contradicted            <NA>  bits=-57.8
#>  outcome_binary            y unverifiable identifiability    unique=2
#>                                                                                  reason
#>            Every observed value is distinct and the column has an identifier signature.
#>  The data cannot say whether a column is the grouping; the entropy balance is reported.
#>  The two columns share no information beyond the independence noise floor (bits=-57.8).
#>                             The data cannot say which two-valued column is the outcome.
#> Structure not in the claim:
#>   paired_pairs       sbp_base, sbp_end (bits=240.1)
#>   coupled_pair       sbp_end, y (pct=69)
```

* **contradicted**: a falsification backed by a value statistic. Here `age` and
  `sbp_end` share no information, so they cannot be a pre/post pair. The exact
  conditions: `subject_id` values repeat while a detected long panel has another
  entity cluster; a paired/agreement pair shares no information, differs more than
  it agrees, or has a constant member; a `coupled_pair` whose members are
  independent; an event column that is not a two-level or 0/1/2 coding, or negative
  times; `repeated_measures` that are not a detected panel's measures, or a wide set
  with no common factor; ordinal `scale_items` with no common factor, or binary items
  forming a one-hot partition; a `time_index` that is not a calendar column; an
  all-distinct or calendar `group_var`; a non-numeric `outcome_continuous`; an
  `outcome_binary` without exactly two observed values. A false alarm costs a
  reviewer one read of the reason; a missed error costs far more.
* **consistent**: the data support the reading, never "correct". For a pair it means
  the two columns share information and their difference is smaller than their
  mean. Binary scale items are never consistent.
* **unverifiable (identifiability)**: design roles (`group_var`, `outcome_*`,
  `covariate`) are not functions of the joint distribution (Lemma 2), so the data
  cannot confirm them.
* **unverifiable (coverage)**: the package has no test for the situation, or the
  evidence is empty (all-NA, fewer than 10 complete cases, non-numeric). No
  time-event association is also coverage, not a contradiction: constant-hazard
  censoring makes the event independent of the observed time.
* Each claimed pair is judged on its own. A survival claim (time plus event) is one
  claim shown as two rows and counted once.
* **Structure not in the claim** lists what the blind detector found that the claim
  does not state under the same role (covariates excluded): here the real pair
  `sbp_base`/`sbp_end`.

## Aligning roles across sites without a codebook

Pooling data from several centers starts with finding which column plays which role
in each file. `harmonize_roles(sites)` runs the name-blind detector on every site and
tabulates, per role, the column each site holds (`aligned` = found everywhere,
`partial`, `absent`). Here the two sites share no column names and the order differs,
and the event is coded 0/1 at one site and 1/2 at the other.

```r
site <- function(seed, nm, n = 200, code = 0:1) {
  set.seed(seed); b <- rnorm(n, 150, 15); t <- rweibull(n, 1.3, exp(5 - rnorm(n)))
  cens <- runif(n, 0, quantile(t, .9))
  d <- data.frame(sample.int(99999, n), sample(c("a", "b"), n, TRUE), b, b - 5 + rnorm(n, 0, 6),
                  round(pmin(t, cens), 1), code[1 + (t <= cens)])
  setNames(d, nm)
}
north <- site(1, c("id", "arm", "sbp_base", "sbp_end", "time", "dead"))
south <- site(2, c("pid", "group", "bp0", "bp1", "ftime", "status"), code = 1:2)[, c(3, 1, 6, 2, 5, 4)]
harmonize_roles(list(north = north, south = south))
#> <role_harmonization> 2 sites: 5 aligned, 0 partial roles
#> 
#>            role  status             north    south
#>       group_var aligned               arm    group
#>    paired_pairs aligned sbp_base, sbp_end bp0, bp1
#>   time_variable aligned              time    ftime
#>  event_variable aligned              dead   status
#>      subject_id aligned                id      pid
#> 
#> Flagged and overlap checks
#>           role         check north south flag
#>      group_var level_jaccard  1.00  1.00 <NA>
#>  time_variable range_overlap  0.72  0.72 <NA>
#> 
#> Descriptive checks (levels, Jaccard, ranges) are in $checks; flags mark structural mismatches only.
```

* `map` is the table above; `checks` adds descriptive rows per role found in two or
  more sites (value type, event levels and coding, group-label Jaccard, range and
  range overlap, battery size); `site_only` lists roles found at one site only.
* Only categorical facts are flagged: `type_mismatch`, `coding_mismatch` (the
  normalized event coding differs: 0/1, 1/2 and TRUE/FALSE normalize to the same
  two-level shape, which level is the event is not checked; a 0/1/2 status differs) and `size_mismatch` (item counts differ). Jaccard and range overlap are
  reported, never flagged, and no numeric threshold decides a flag.
* Alignment is by role only; names appear in the output cells. The residual
  `covariate` role is excluded.
* `aligned` means the detector found the role at every site; it does not establish that the columns measure the same thing; read `checks` (levels, Jaccard, ranges) before pooling. Cell order is alphabetical and carries no pairing direction.

## Known limitations

- A **numeric-coded** treatment (0/1) is not distinguishable from a numeric covariate
  from the data alone; `group_var` is detected for character/factor groupings.
- A binary outcome with a strictly positive continuous cause can carry the same value
  signature as a survival time/event pair; the package may then report a survival pair.
- When a positive, right-skewed column is strongly associated with a two-valued column
  that is not the documented event (e.g. `survival::veteran`, `diagtime`/`prior`), that
  pair is reported instead of the documented one. Name-blind, the two cannot be told apart.
- The binary-outcome evidence is a two-part MDL code (0.3.0): each admitted cause
  must pay its parameter bits plus the bits naming it among the m candidates, which
  keeps phantom outcomes on pure noise to roughly 6% of tables -- at the price that
  a weak, borderline true cause may stay honestly silent.
- When the paired-measurement scorer latches onto a (cause, outcome) column pair,
  the continuous outcome can be claimed by `paired_pairs` and stay silent; lifting
  that exclusion was measured to misassign more outcomes than it recovered.
- Paired versus agreement is the sign of the MDL shift bits of the single best pair, a
  BIC-type test of a mean difference between the two columns. `agreement_pairs` therefore
  means "no detectable systematic shift": biased method-comparison (Bland-Altman) data with
  a real calibration bias is labeled `paired_pairs`, and the label is n-dependent (a
  0.5-unit bias at SD 15 flips from agreement to paired between n of about 200 and 1000).
  About 3% of unbiased pairs at n = 100 read as paired under the null. Exactly one pair is
  reported across the two roles.
- `name_bonus` adds 10 points on the 100 * B group scale, where candidates usually differ
  by 1-5 points, so a group keyword hint is decisive in practice (a no-op under `NULL`).
- Unequal allocation: a 2:1 arm can lose the group ranking to a more balanced binary
  covariate (a documented ceiling of the entropy ranking).
- Survival, wide `repeated_measures` and binary-item `scale_items` still use the ported
  point scorers (`pct >= 30`). Survival: ranking pairs by association bits picks the
  strongest dependence (on `survival::pbc`, bilirubin/status), and every variant that
  repaired this added false positives to the survival roles. Wide repeated measures: a
  one-factor code on standardized columns measures shared variance only and discards the
  same-units structure (equal means and spreads across repeats) that the scorer uses; four
  block-growth variants raised false positives from 15 to 163-207 tables (dev-band measurement, pass 5). Binary items: an
  all-pairs dependence test rejected every simulated binary battery.
- For ordinal batteries of k <= 3 items the one-factor model is saturated (k = 2) or just-identified (k = 3), so the
  test there is only whether the items co-vary; unidimensionality is enforced for k >= 4.
- The ordinal `scale_items` candidate set is not pruned: a stray ordinal covariate can ride
  along with a real battery.
- A categorical column in which any level is observed fewer than 3 times is excluded from
  `group_var` candidacy.
- Accuracy has been confirmed on a synthetic data-generating process only; a
  pre-registered evaluation on held-out real data is planned.

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
