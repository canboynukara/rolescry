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

## Citation

* rolescry now has its own archival DOI, 10.5281/zenodo.21003941;
  `citation("rolescry")` and `CITATION.cff` reference both rolescry and the
  parent MDStatR engine.
