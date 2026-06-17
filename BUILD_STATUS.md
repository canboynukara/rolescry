# rolescry — Build Status (Phase B)

**Package:** `rolescry` 0.1.0 — Name-Blind Variable-Role Detection by Data Signature
**Built:** 2026-06-17 · **R:** 4.5.3 (ucrt), Windows 11 · **From:** MDStatR detector subsystem (copy-only)
**Location:** `C:\Users\cboyn\Desktop\rolescry\` · **Tarball:** `rolescry_0.1.0.tar.gz`

---

## 1. `R CMD check --as-cran` verdict

```
Status: 1 NOTE, 0 ERRORs, 0 WARNINGs   (exit 0)
```

Test suite: **all pass** (`checking tests ... OK`; 8 test files, ~107 assertions, run on the
*installed* tarball inside the check). Vignette builds and re-builds OK. Examples run OK.
Source is ASCII-clean. S3 generic/method consistency OK. No unstated dependencies.

### The single residual NOTE (pre-publication, expected)

```
checking CRAN incoming feasibility ... NOTE
Found the following (possibly) invalid URLs:
  URL: https://github.com/canboynukara/rolescry           Status: 404
  URL: https://github.com/canboynukara/rolescry/issues    Status: 404
```

**Rationale:** the GitHub repository does not exist *yet*. These are the intended canonical
`URL`/`BugReports` for the package; the 404 resolves the moment the operator creates the repo
(publication step 1 below). This is the only residual and it is not a code defect.

### Environment notes (cleared, for the record)

- **LaTeX absent locally** → checked with `--no-manual` (the PDF reference manual; the `.Rd`
  files themselves are fully checked and pass). CRAN/r-universe build the manual.
- **pandoc** is not on `PATH` by default on this machine; it was supplied via
  `RSTUDIO_PANDOC` + prepended to `PATH` (RStudio's bundled Quarto pandoc) so the vignette
  builds and the README/NEWS check passes. CRAN/r-universe have pandoc, so this is automatic
  there.

---

## 2. Severed couplings (9 bonds → all cut)

| # | MDStatR bond | How it was severed in `rolescry` |
|---|---|---|
| 1 | `%||%` injected global | Defined once, internal (`R/utils.R`), NA-aware variant. |
| 2 | `MDSTATR_VERSION` global + load banner | Load banner deleted; no global referenced. |
| 3 | `md_casefold_ascii` (stringi + `Sys.setlocale`) | Replaced by base-only `.casefold_ascii` (`chartr`, no locale mutation). |
| 4 | `classify_value_type` via `source("engine/...")` ×3 | Bundled into `R/classify.R`; all 3 `source()`/`file.exists()` lines removed. |
| 5 | `score_id_candidate` `exists()`-probe | Not ported (optional in source); value-based `.build_var_info` used instead. |
| 6 | Ax95 harness wrapper (`run_097_detect_roles`) | Dropped; pure `run_all_detections` exposed via `detect_roles()`. |
| 7 | `reason_code_registry.csv` coupling | Gone with the Ax95 wrapper; package returns a plain S3 object. |
| 8 | `relationships` arg (never built upstream) | Not needed; the NMI/entropy math kept, relationship penalties dropped. |
| 9 | globalenv source-once sentinels (P-Audit-5) | Stripped entirely; package collation loads each file once. |

Additional severance: the Shiny `rv` object was replaced by an explicit `var_info` data.frame
parameter (the Phase A baseline ran with `config = NULL`/`confirmed_roles = NULL`, entering the
mathematical layer directly, so this is exactly equivalent).

## 3. CRAN-blocker resolutions (180 Phase A findings → resolved)

| Category | Phase A count | Resolution |
|---|---:|---|
| `unqualified_call` | 99 | Every base-extra call fully qualified (`stats::`, `utils::`). |
| `non_ascii` | 33 | All comments/strings translated to ASCII English; `\p{L}` (perl) replaces the literal Latin/Cyrillic letter range. Source is ASCII-clean (`checking ... non-ASCII characters ... OK`). |
| `console_write` (`cat`/`print`) | 18 | Replaced by `.say()` (verbose-gated `message()`); load banners deleted. |
| `locale_dictionary` | 13 | Externalized into the optional `name_bonus` (default `NULL`), ASCII-English, via `rolescry_default_name_bonus()`. |
| `global_state` | 9 | `%||%` internalized; `MDSTATR_VERSION` dropped; globalenv sentinels removed; `exists()` probes removed. |
| `filesystem` (`source`/`inst/*.json`) | 5 | All removed; package uses NAMESPACE collation, not runtime `source()`. |
| `other` (load `message`, `Sys.setlocale`) | 3 | Load messages removed; locale-mutating casefold replaced. |
| `internet` / `TF_literal` / `seq_1_to_length` | 0 | Already clean; preserved. |

## 4. Baseline-equivalence & H4 turnusol invariant

- **H4 (keystone):** a `col_N` dataset and its named twin (byte-identical data) yield
  **identical role assignments by column index**. Verified in `tests/testthat/test-turnusol.R`
  and re-confirmed on the installed package (`H4 ... TRUE`).
- **Baseline-equivalence:** the signature detectors reproduce the Phase A Pipeline-1 behavior
  exactly. Confirmed during development against MDStatR DS27/DS40/DS41
  (`DS27_equiv=TRUE, DS40_equiv=TRUE, H4_invariant=TRUE`), and pinned in
  `tests/testthat/test-baseline-equiv.R` on seeded in-test fixtures (no MDStatR CSVs copied).
- **Name-bonus path** (the Phase A "unexercised" caveat) is now exercised:
  `tests/testthat/test-name-blind.R` shows `name_bonus` tie-breaks a binary outcome to `death`
  on the named twin while remaining inert on the `col_N` twin.
- **Determinism:** no RNG in the detector; identical output across runs.

## 5. Package contents

```
rolescry/
  DESCRIPTION  NAMESPACE  LICENSE.md  README.md  CITATION.cff  .Rbuildignore
  R/      rolescry-package.R utils.R classify.R var_info.R signatures.R nmi.R
          detect.R detect_roles.R name_bonus.R header.R read_data.R
  man/    detect_roles.Rd detect_header.Rd read_data.Rd compute_nmi.Rd
          rolescry_default_name_bonus.Rd rolescry-package.Rd
  tests/testthat/  test-turnusol.R test-baseline-equiv.R test-name-blind.R
          test-classify.R test-nmi.R test-header.R test-detect-roles.R
          test-read-data.R helper-dgp.R
  vignettes/rolescry.Rmd
  inst/CITATION
```

**Public API (exported):** `detect_roles()`, `read_data()`, `detect_header()`, `compute_nmi()`,
`rolescry_default_name_bonus()`, plus S3 `print`/`summary` for `role_detection`. Everything
else is internal.

## 6. Citation & attribution

- `CITATION.cff` (root), `inst/CITATION` (R-native), and README all credit the parent work.
- README attribution line: *Derived from Boynukara, C. (2026). MDStatR (v2.1.0 Veritas).
  Zenodo. https://doi.org/10.5281/zenodo.20707791*.
- `citation("rolescry")` returns two entries (the package + the parent MDStatR Zenodo record).
- License: Apache-2.0 (inherited from MDStatR); `LICENSE.md` ships the full text.

## 7. Remaining manual steps to publish (operator-only — no git/remote done here)

1. **Create the GitHub repo** `https://github.com/canboynukara/rolescry`, then in the package dir:
   `git init && git add -A && git commit -m "rolescry 0.1.0" && git branch -M main`
   `git remote add origin https://github.com/canboynukara/rolescry.git && git push -u origin main`
   (This resolves the only residual `R CMD check` NOTE — the URL 404.)
2. **Tag the release:** `git tag v0.1.0 && git push --tags` (matches `CITATION.cff` version 0.1.0).
3. **Register on r-universe:** add `{ "package": "rolescry", "url": "https://github.com/canboynukara/rolescry" }`
   to your `<user>.r-universe.dev` registry (the `packages.json` in your r-universe `universe`
   repo). r-universe then builds binaries; users install with
   `install.packages("rolescry", repos = "https://canboynukara.r-universe.dev")`.
4. *(Optional, for Zenodo archiving of the new package)* enable the GitHub–Zenodo integration and
   publish a release to mint a DOI for `rolescry` itself.

No `git`, push, or remote action was performed by this build — publication is left to the operator.

## 8. Source-repo safety attestation

`C:\Users\cboyn\Desktop\MDStatR` was accessed **read-only** (Read/Grep/Glob + `Copy-Item` *from*
it). No writes, edits, moves, deletes, or git operations were performed against it. All 10 source
and dataset files inspected retain their pre-session modification times (0 modified during this
session). All generated artifacts live under `C:\Users\cboyn\Desktop\rolescry\` (and the Phase A
report/sandbox under `C:\Users\cboyn\Desktop\PHASE_A_REPORT\`).
