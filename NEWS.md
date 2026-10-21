# rolescry 0.3.0

Corrections so that the package does what its documentation says, found by running it on
public survival data. Detection stays deterministic and name-blind.

## New function: `verify_roles()`

* `verify_roles(data, claim)` checks a declared assignment of roles against the
  name-blind value signature. Verdicts: `contradicted` (a falsification backed by a
  value statistic), `consistent` (a reading the data support, never "correct") and
  `unverifiable` with sub-reason `identifiability` (design roles: group, outcome,
  covariate) or `coverage` (no test exists, or the evidence is empty). Binary scale
  items are never consistent. `contradicted` is issued only for falsifications backed
  by a value statistic (exact list in `?verify_roles`); absence of a time-event
  association is `unverifiable(coverage)`, not a contradiction. Each claimed pair is
  judged on its own, and a survival claim (time plus event) is one claim shown as two
  rows. `unclaimed` lists blind-found structure the claim does not state under the
  same role.
  Detectors are unchanged. See `?verify_roles`, the README and the vignette.

## New function: `harmonize_roles()`

* `harmonize_roles(sites)` aligns the roles detected name-blind in two or more sites
  (a named list of data.frames) into a `map` (role by site, `aligned`/`partial`/`absent`),
  descriptive `checks` and a `site_only` table. Flags are categorical only:
  `type_mismatch`, `coding_mismatch`, `size_mismatch`; Jaccard and range overlap are
  reported but never flagged. The residual `covariate` role is excluded. Detectors
  are unchanged. See `?harmonize_roles`, the README and the vignette.

## New role: `coupled_pair`

* The identifiable object in a dependent (binary, positive-continuous) pair is the
  PAIR itself; whether it is (survival time, event) or (cause, binary outcome) is
  design-relative and, in general, not a function of the joint distribution. The
  survival roles are now assigned only when censoring geometry singles the survival
  reading out: the continuous member is anchored at zero (mean - 2 sd <= 0), events
  sit strictly earlier than censorings, and the top time-tertile's event rate is
  strictly below the bottom's. Otherwise the new `coupled_pair` role reports the
  pair with both readings, and the user decides. Documented ceiling: a zero-anchored
  continuous cause whose larger values are protective still reads as survival.
* A randomized-trial table (binary arm shifting a positive measurement) is therefore
  no longer reported as a survival pair: it becomes `coupled_pair`.
* **Dominance gate.** When two decisive pairs on DISJOINT time columns compete,
  the winner is assigned the survival roles only if its event-time association
  dominates at the 2-sigma universal (Fisher z two-correlation comparison);
  otherwise the pair is reported as `coupled_pair` with the competing reading
  named. A rival sharing the winner's time column agrees on the axis and merely
  offers another event candidate (a ranking question), and an ID-typed time
  candidate may win on its exceptional evidence bar but may not veto a
  structurally-typed reading.

## Binary-outcome evidence: two-part MDL code

* The per-predictor 2-sigma dCor test behind `outcome_binary` is replaced by a
  Bernoulli MDL code: a candidate cause is admitted only when its one-slope logistic
  code -- or, for non-monotone dependence, a two-Gaussian split of the predictor by
  the outcome -- saves more bits than the slope's parameter cost (1/2 log2 n) PLUS
  the log2(m) model-index bits naming it among the m candidates (Rissanen two-part
  code). The naming term makes the evidence selection-consistent, so scanning many
  predictors on noise no longer buys phantom outcomes (~14% -> ~6% of pure-noise
  tables); a weak borderline cause now stays honestly silent instead.

## Pass 5: derived deciders for group, pairs and ordinal scales

* **`group_var`** is the character/factor column with the highest normalized Shannon
  entropy H/log2(k) (ties: fewer levels, then value fingerprint); a column with a level
  seen fewer than 3 times is not a candidate. The point table is gone.
* **Paired vs agreement** is decided by the MDL shift bits of the best pair (bits saved by
  coding the difference with a free mean, net of one parameter): > 0 => `paired_pairs`,
  otherwise `agreement_pairs`. `potential_pairs` now carries `bits` and `shift_bits`.
* **Ordinal `scale_items`** are accepted when a one-factor Gaussian two-part code on the
  standardized candidate columns saves bits over independence (`.mdl_factor_bits`).
* `score`/`max_score`/`pct` are bits (`max_score` NA) for derived roles and points for the
  ported ones.
* **Measured and rejected, kept on ported scorers:** survival admission/rank by association
  bits, wide `repeated_measures` as a one-factor block (false positives rose from 15 to
  163-207 tables), and a derived binary-item scale criterion (it rejected every simulated
  binary battery) each regressed false positives or recall and were reverted. Survival, wide `repeated_measures` and binary `scale_items` therefore
  remain on point scorers.

## Domain generalization (finance and beyond)

* **New role `time_index`.** Date/POSIXct columns -- and character columns whose
  every value parses as a date -- type as Datetime and are reported as calendar
  axes. They were previously typed "Other" and invisible to every detector. The
  class and format are properties of the values, so name-blindness is preserved.
* **Panel entity/time axis.** On a balanced long panel both axis readings satisfy
  the structural test. Serial dependence of the measures along the candidate time
  axis (2-sigma) now picks the entity; when nothing discriminates, the panel is
  reported with its admissible axis readings (`repeated_measures$panel_axes`) and
  `subject_id` stays silent. A Datetime column can be the panel time axis.
* **Short rating batteries.** A battery of exactly 3-4 ordinal items (e.g.
  multi-agency credit ratings) is now scoreable as `scale_items`, guarded by the
  pairwise 2-sigma dependence gate; the >= 5-item path is unchanged.
* **Finance keyword dictionary.** `rolescry_default_name_bonus(domain = "finance")`
  (or `"all"`); the default `"clinical"` is unchanged. The name channel remains
  optional and capped.
* README, DESCRIPTION and examples reworded domain-neutrally.

## Behaviour changes

* **Logical event columns.** A TRUE/FALSE column now value-types as Binary; it was
  typed "Other" and silently excluded from every survival and outcome candidate set.
* **Sink veto in repeated measures.** A column MDL-explained by two or more mutually
  independent other block members is a dependency sink (collider), not a repeated
  measurement, and is removed from the block. `y = x1 + x2 + noise` is no longer
  swallowed as a correlated block, so the continuous outcome can be detected.

* **`covariate` is now the residual role.** It lists the continuous, categorical and
  binary columns that no other role claimed. In 0.2.0 it listed every such column,
  including the survival pair, the group and the outcomes.
* **Survival pairs are compared on common rows.** Every admissible time/event pair is
  re-evaluated on the rows complete in all candidate columns. Previously each pair was
  judged on its own available cases, so a pair whose columns were missing in different
  rows could win on a different subsample (e.g. `survival::lung`).
* **Integer follow-up times.** An all-distinct positive integer column (which value-types
  as an identifier) is now a survival-time candidate. Because this overturns the
  identifier typing, it must pass both association tests (monotone MDL and the dCor
  independence test), not either one (e.g. `survival::ovarian`).
* **Subject id atomicity.** A column already assigned to a survival role is no longer
  also reported as `subject_id`.

## Bug fixes

* The KSG mutual-information estimator returned `NaN` on tied values (a neighbour count
  of -1 passed to `digamma()`), which made the dependence measure `NaN` and emitted
  warnings (e.g. `survival::pbc`). Counts are now clamped at zero.
* Correlation calls on constant columns in the paired/agreement scorers no longer emit
  "standard deviation is zero" warnings; the result is unchanged (`NA`).

## Documentation

* The README quick-start example now uses data in which the documented roles are
  actually present, and a "Known limitations" section was added.
* `?detect_roles` states that roles can overlap, that `covariate` is residual, and that
  a role is reported as not found rather than guessed.

# rolescry 0.2.0

Information-theoretic re-founding of the detection engine, a correctness fix to
the name-blindness contract, and more robust survival-code handling. Detection
remains deterministic and, by default (`name_bonus = NULL`), a pure function of
each column's value signature.

## Behaviour changes

* **Name-blindness contract fix.** Column typing no longer consults column names
  when `name_bonus = NULL`. Previously a weak-signature identifier column could
  be typed from its name; such a column is now typed from its values only, so
  its assignment can differ from 0.1.0. This restores the package's core
  invariant: renaming columns to `col_1, col_2, ...` cannot change any role.
* **Survival event coding.** Event detection now normalises 2- and 3-level
  numeric and logical codings to `{0, 1}`; factor and character codings are
  intentionally not coerced.

## Detection improvements

* Dependency structure is scored with normalized mutual information and
  distance-based criteria, so nonlinear dependence is visible to the detector.
* Continuous and binary outcome detection is recast as a dependency *sink*: a
  column explained by two or more mutually independent predictors. The detector
  stays silent when no column is decisively a sink, rather than guessing.
* Repeated-measures detection is set-based, and every single-winner selection
  uses value-fingerprint canonical tie-breaks, so results do not depend on
  column order.

## Invariance guarantees

* RELABEL (renaming every column) and S_n (permuting columns) both leave every
  role assignment unchanged (|delta| = 0). These are structural properties of
  the estimator and hold independently of any tuning.

## Validation and scope

* A single pre-registered, held-out confirmatory run (OSF registration
  osf.io/8ecau; analytic-code archive DOI 10.5281/zenodo.21003941) met its
  binding criteria on the specified synthetic data-generating process:
  structural invariance (RELABEL and S_n |delta| = 0), admissible macro-F1 above
  the pre-registered floor, and superiority over the strongest baseline.
* Scope: this confirms estimator stability on that synthetic process only. It
  does **not** establish real-data or external validity, which is a separate,
  separately registered question. Two limitations are disclosed: `event_variable`
  fell just outside its pre-registered calibration band, and `covariate` remains
  inadmissible (high false-positive rate) and is excluded from the headline.
* Packaging only (no algorithm change): the CRAN 0.2.0 source differs from the
  archived analytic snapshot (git tag `v.0.2.0`; Zenodo 10.5281/zenodo.21003941)
  in corrected vignette example data, test expectations aligned with the current
  detector, a `stats::` import, and an updated `.Rbuildignore`. Detection
  behaviour, its thresholds, and the confirmatory results above are unchanged.

## Citation

* rolescry now has its own archival DOI, 10.5281/zenodo.21003941;
  `citation("rolescry")` and `CITATION.cff` reference both rolescry and the
  parent MDStatR engine.
