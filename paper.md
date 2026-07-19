---
title: "rolescry: Name-blind detection of statistical variable roles by information-theoretic data signature in R"
tags:
  - R
  - exploratory data analysis
  - information theory
  - minimum description length
  - reproducibility
  - variable role detection
authors:
  - name: "Can Boynukara, MD"
    orcid: 0000-0002-6075-3923
    affiliation: 1
affiliations:
  - name: "Department of Internal Medicine, Acıbadem University, Istanbul, Türkiye"
    index: 1
date: 19 July 2026
bibliography: paper.bib
---

# Summary

`rolescry` is an R package that assigns statistical *roles* to the columns of a
data frame from the information-theoretic signature of their values and their
inter-column dependency structure, rather than from their names. Given a data
frame, it labels columns with one or more of eleven roles — group variable,
paired and agreement measurements, survival time and event, subject identifier,
continuous and binary outcome, repeated measures, scale items, and covariate —
and remains silent where a role is not identifiable from the data alone. The
governing principle is *Data inspice, non nomen* ("inspect the data, not the
name"): renaming every column to `col_1, …, col_N` leaves every assignment
unchanged (relabelling invariance), and permuting the column order permutes the
output identically (permutation equivariance). Detection is fully deterministic
and runs offline: it uses no large language models and transmits no data
externally. Column names influence the result only through an optional, capped,
off-by-default hint channel.

# Statement of need

Analysts routinely need to know *what statistical part each column plays* before
any modelling begins — which column is the grouping factor, which is the
outcome, which pair encodes survival time and event, which columns are repeated
measures of one subject. In R this step is currently either manual or
name-driven. Profiling tools such as `skimr` [@skimr], `DataExplorer`
[@DataExplorer], and `dataMaid` [@dataMaid] summarise and type-check columns but
do not infer their statistical role; modelling frameworks such as `recipes`
[@recipes] require the *user* to declare roles. The closest automated work,
semantic-type detection (`Sherlock` [@hulsebos2019], `Sato` [@zhang2020sato]),
targets real-world *entity types* (e.g. name, address) using machine learning on
column names and surrounding context in Python — it is name- and
context-dependent and does not recover *statistical* roles.

`rolescry` fills this gap with a deterministic, name-blind, information-theoretic
inference of statistical roles. Its distinguishing property is a structural
guarantee rather than a benchmark score: because every estimator is a function
of permutation- and relabelling-invariant column-value features, the assignment
is invariant to column renaming and equivariant to column reordering by
construction. This makes the tool useful for pre-analysis screening, for
scaffolding reproducible analysis pipelines, for generating name-independent
dataset metadata, and for teaching, in settings where reliance on column names
is fragile or where offline, deterministic, privacy-preserving operation is
required.

# How it works

`rolescry` reads each column's information-theoretic signature [@shannon1948] and
the dependency structure between columns in four layers, so that *what* it
detects is tied to *how* the data compresses rather than to any label:

1. **Dependency backbone — what depends on what.** Pairwise dependence is
   measured by distance correlation [@szekely2007], which is zero if and only if
   two columns are independent and, unlike linear correlation, captures nonlinear
   and non-monotone relations. This builds a dependency graph over columns from
   value distances alone, so it is unchanged by renaming or reordering.

2. **Role scoring by compression — why a column earns a role.** A candidate role
   is asserted only if its model *compresses* the column below the null model,
   using a two-part minimum description length criterion [@rissanen1978;
   @grunwald2007]. Decision thresholds are therefore *derived* from sample size
   (the `(k/2)·log₂(n)` parameter cost), not hand-tuned — collapsing seven tuned
   constants of an earlier points-based system to zero while preserving the
   literature-based small-sample floors. Concretely: an outcome is a dependency
   *sink* explained by two or more mutually independent predictors; a subject
   identifier is a near-unique integer signature; scale items are a correlated
   block sharing one latent factor; a survival pair is a time column whose events
   compress it.

3. **Direction, only where identifiable.** For a detected continuous outcome,
   cause→effect direction is annotated by additive-noise-model asymmetry
   [@hoyer2009; @peters2014]; in the linear-Gaussian symmetric case direction is
   unidentifiable, so `rolescry` stays *silent* rather than guessing.

4. **Confidence and silence.** Every assignment carries a description-length
   "bit-margin" — how many bits separate it from the next-best explanation — and
   the tool abstains when that margin is within coding noise.

Ties between equally good columns are resolved by an order-free value
fingerprint; this is what makes the output invariant to renaming (relabelling
Δ = 0) and equivariant to reordering (permutation Δ = 0). Both properties hold by
construction and are re-certified empirically at 20 random relabellings and
permutations per dataset. A worked, layer-by-layer walk-through is provided in
the package vignette.

# Validation

Estimator behaviour was assessed by a single pre-registered confirmatory run
[@nosek2018] on a frozen synthetic data-generating process (OSF registration
`osf.io/8ecau`; analytic code archived at Zenodo [@rolescryZenodo]). On 200
held-out synthetic datasets the package attained an admissible macro-F1 of 0.67
(95% CI [0.63, 0.70]), above its pre-registered floor of 0.61 and the best
deterministic baseline (0.22), with the structural invariance criterion passing
on every dataset. **Scope (stated deliberately):** this confirms the estimators'
behaviour on the specified synthetic data-generating process; it does *not*
establish real-data or external validity, which is out of scope here and is the
subject of a separate, separately pre-registered study.

# Availability

`rolescry` is distributed on CRAN under the Apache-2.0 licence, with source and
issue tracker in the project repository and an archived release at Zenodo
[@rolescryZenodo]. The package depends only on base R and its recommended
packages.

# Acknowledgements

`rolescry` is developed by a single author; benchmarking was author-conducted and
its limitations (author-designed oracle, synthetic-only confirmatory scope) are
documented in the validation materials. No funding is reported.

# References
