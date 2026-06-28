# Role detection engine. Ported from MDStatR engine/run_097_detect_roles.R
# (detect_variable_role + .detect_* + run_all_detections), with the Shiny rv
# coupling severed: var_info is an explicit data.frame argument, and the
# confirmed_roles / config layers (Shiny-runtime concerns) are dropped -- the
# Phase A baseline ran with config = NULL and confirmed_roles = NULL, so it
# already entered the mathematical layer directly, making this exactly
# equivalent. The optional capped name bonus is threaded into the
# candidate-selecting detectors (group, outcome, subject id, survival); with
# name_bonus = NULL it is a no-op (purely mathematical). No exports.

.detect_group_var <- function(data, var_info, name_bonus = NULL) {
  # D-A (council ruling, Tier-3 boundary): group_var conceptually SPLITS into (i) a DESIGNED,
  # balanced-exogenous RCT arm -- balanced vs EVERY baseline covariate (group _|_ covariates,
  # the do-operator's severed incoming edge) AND directed-upstream of the outcome (group _|/_
  # outcome) -- identifiable only in construction-true RCT data; vs (ii) an OBSERVATIONAL /
  # numeric-coded treatment, which is estimand-relative (confounded with covariates) and is
  # EXCLUDED from accuracy, routed to honest SILENCE exactly like `confounder` -- never a
  # fabricated pick. Detector (i) is GATED behind a generator-D RCT-arm extension and is NOT
  # enabled here until that construction-true signal exists and fires ONLY on the RCT arm
  # (zero firing on a balanced-observational arm -- else it would breach the Gödel/Hume
  # firewall). Until then group_var stays char/factor-only (the twice-reverted numeric wall).
  # C3 finding (REVERTED): extending candidates to low-cardinality NUMERIC columns to
  # catch numeric-coded treatment arms regressed sim group_var 0.672 -> 0.422 (precision
  # 0.33) with NO real gain -- numeric low-card columns (scale items, balanced binary
  # outcomes/events) flood in as false positives, and a conditioning/source guard cannot
  # separate them (they all condition something). Name-blind, a numeric-coded treatment
  # is not distinguishable from a low-card outcome/ordinal/scale-item without external
  # design information -- a Tier-3 identifiability limit. Kept char/factor-only.
  candidates <- var_info$column[var_info$user_type %in% c("Categorical", "Binary")]
  cols <- character(0); ranks <- numeric(0); gss <- list()
  for (col in candidates) {
    gs <- .score_group_signature(data[[col]])
    if (gs$pct >= 30) {                                  # LEGIT-threshold preserved
      cols <- c(cols, col)
      ranks <- c(ranks, gs$pct + .name_bonus_pts(col, name_bonus$group_var))
      gss[[length(gss) + 1]] <- gs
    }
  }
  # D3 disciplined re-attempt (REVERTED): the SOURCE-vs-SINK directional asymmetry --
  # admit a low-card NUMERIC column iff it is exogenous (NOT MDL-compressed by the others,
  # so not an outcome SINK) AND high fan-out (.dep-conditions >= 2 columns) -- was built
  # and MEASURED. It regressed sim group_var 0.6718 -> 0.6255 (precision 0.567) AND real
  # 0.0667 -> 0.0000 (dangerous-class). Mechanism: the asymmetry separates a group from an
  # OUTCOME (a sink) but NOT from a numeric COVARIATE (both are exogenous high-fan-out
  # SOURCES), so used covariates flood as false positives. Reverted per the firewall; this
  # CONFIRMS (a 2nd time, after the C-line) the name-blind group/covariate identifiability
  # wall: a numeric-coded treatment is indistinguishable from a low-card covariate without
  # external design info (Tier-3-adjacent). group_var stays char/factor-only.
  pick <- if (length(cols) > 0) .canonical_pick(cols, ranks, data) else NA_character_
  if (!is.na(pick)) {                                    # NA = unbreakable tie -> silence (C1)
    gs <- gss[[match(pick, cols)]]
    return(list(
      found = TRUE, role = "group_var", column = pick,
      detected_by = if (is.null(name_bonus$group_var)) "mathematical" else "mathematical+name",
      score = gs$score, max_score = gs$max, pct = gs$pct, breakdown = gs$breakdown,
      details = gs$details
    ))
  }
  list(
    found = FALSE, role = "group_var", column = NULL, detected_by = "mathematical",
    score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

.detect_pairs_generic <- function(data, var_info, scorer, role) {
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  if (length(num_cols) < 2) {
    return(list(found = FALSE, role = role, column = NULL, detected_by = "mathematical",
                score = 0, max_score = 110, pct = 0, breakdown = list()))
  }
  pairs <- list()
  for (i in seq_along(num_cols)[-length(num_cols)]) {
    for (j in (i + 1):length(num_cols)) {
      a <- data[[num_cols[i]]]; b <- data[[num_cols[j]]]
      # C2: a measurement pair is decided by the DERIVED MDL criterion -- the
      # consensus+difference factorization compresses (bits > 0 => the columns share
      # information) AND the STRUCTURAL agreement condition var(a-b) < var((a+b)/2)
      # ("agree more than they differ"). This REPLACES the pct>=30 accept, the
      # paired/agreement r-windows, and the dCor 0.30 gate with no tuned constant. The
      # .mdl_pair_bits test is O(n), so it runs on all pairs cheaply; the scorer (kept
      # only for a descriptive, NON-deciding breakdown) is computed on survivors.
      mb <- .mdl_pair_bits(a, b)
      # D1 Gaussian-reduction guard. The pair bits = (n/2)log2(va*vb/(vm*vd)) equal
      # n*MI(a,b) EXACTLY only for a Gaussian-admissible pair (the measurement-pair
      # model: two noisy readings of one latent -> bivariate Gaussian). The closed-form
      # acceptance is mb$bits>0 && mb$agree. Where both margins pass `.gauss_margin`,
      # that closed form IS the exact NMI test. Where a margin is non-Gaussian the
      # variance bits can UNDERSTATE the true MI, so the KSG/dCor estimator .dep is
      # added as a fallback that can only ADD a non-Gaussian agreeing pair the closed
      # form rejected (clearing the 2/sqrt(m) independence floor -- the SAME 2-sigma
      # universal as the ANM margin, no new tuned constant). Ranking ALWAYS uses the
      # variance-MDL bits mb$bits, so every pair the closed form accepts is IDENTICAL to
      # the frozen behaviour (a NO-OP on the linear-Gaussian sim); the fallback only
      # appends off-Gaussian catches (mb$bits<=0, ranked last). The STRUCTURAL pair
      # condition var(d) < var(m) (mb$agree) gates both routes.
      cf <- isTRUE(mb$bits > 0 && mb$agree)
      accept <- cf
      if (!cf && !(.gauss_margin(a) && .gauss_margin(b))) {
        m_eff <- sum(is.finite(a) & is.finite(b))
        accept <- isTRUE(mb$agree) && .dep(a, b) > 2 / sqrt(max(m_eff, 1))
      }
      if (accept) {
        pairs[[length(pairs) + 1]] <- list(col1 = num_cols[i], col2 = num_cols[j],
                                           score = scorer(a, b), bits = mb$bits)
      }
    }
  }
  if (length(pairs) > 0) {
    # canonical, order-free pick: MDL bits desc, then symmetric value-fingerprint asc
    bts <- vapply(pairs, function(p) p$bits, numeric(1L))
    fps  <- vapply(pairs, function(p) paste(sort(c(.col_fingerprint(data[[p$col1]]),
                                                   .col_fingerprint(data[[p$col2]]))), collapse = "~"),
                   character(1L))
    fps <- unname(fps); ord <- order(-bts, fps)
    tie <- length(ord) >= 2L && isTRUE(bts[ord[1]] == bts[ord[2]]) && isTRUE(fps[ord[1]] == fps[ord[2]])
    if (!tie) {                                            # NA-silence on value-identical pair tie
      top <- pairs[[ord[1]]]
      return(list(
        found = TRUE, role = role,
        column = list(list(col1 = top$col1, col2 = top$col2)),
        all_pairs = pairs[ord], detected_by = "mathematical",
        score = top$score$score, max_score = top$score$max, pct = top$score$pct,
        breakdown = top$score$breakdown
      ))
    }
  }
  list(found = FALSE, role = role, column = NULL, detected_by = "mathematical",
       score = 0, max_score = 110, pct = 0, breakdown = list())
}

.detect_paired_pairs <- function(data, var_info) {
  .detect_pairs_generic(data, var_info, .score_paired_signature, "paired_pairs")
}
.detect_agreement_pairs <- function(data, var_info) {
  .detect_pairs_generic(data, var_info, .score_agreement_signature, "agreement_pairs")
}

.detect_survival_components <- function(data, var_info, component, name_bonus = NULL) {
  binary_cols <- var_info$column[var_info$user_type %in% c("Binary", "Categorical")]
  numeric_cols <- var_info$column[var_info$user_type %in% c("Continuous", "Ordinal")]
  # C3: also admit low-cardinality integer status columns (e.g. pbc {0,1,2},
  # which value-types as Continuous), normalized to {0,1} by .normalize_event.
  int_status <- numeric_cols[vapply(numeric_cols, function(c) {
    v <- suppressWarnings(as.numeric(data[[c]])); v <- v[!is.na(v)]
    length(v) > 0 && all(abs(v - round(v)) < 1e-9) && length(unique(v)) == 3L
  }, logical(1L))]
  event_cands <- unique(c(binary_cols, int_status))
  evs <- character(0); tms <- character(0); ranks <- numeric(0); sss <- list()
  for (ev in event_cands) {
    ev01 <- .normalize_event(data[[ev]])                 # C3: {1,2}/factor/logical/{0,1,2} -> {0,1}
    if (is.null(ev01)) next
    obs <- !is.na(ev01)                                  # C2: available-case, no hard NA skip
    if (sum(obs) < 10) next                              # n<10 LEGIT-threshold preserved
    if (length(unique(ev01[obs])) != 2L) next
    for (tm in numeric_cols) {
      if (tm == ev) next
      tm_vals <- suppressWarnings(as.numeric(data[[tm]]))
      if (any(tm_vals <= 0, na.rm = TRUE)) next
      ss <- .score_survival_signature(tm, ev, data)
      if (ss$pct >= 30) {                              # LEGIT-threshold preserved
        # A GENUINE survival pair has the event statistically associated with the
        # time (events cluster at different times than censoring); a spurious
        # (balanced-binary, positive-continuous) pair does not. Require & reward
        # that association so e.g. status+time beats sex+meal.cal. Value-based =>
        # equivariant.
        cc <- obs & !is.na(tm_vals)
        # Survival association: the linear-MDL break-even floor (event compresses time
        # under a 1-predictor Gaussian model) + |cor| ranking weight. C1 tried a
        # variance-aware 2-group MDL to also catch NON-monotone hazard, but its total
        # bits scale with n, so on large-n monotone REAL survival it over-admitted
        # non-survival binary-vs-continuous variance-splits and regressed event_variable
        # (0.75 -> 0.44). D3 resolves this WITHOUT that regression: the linear MDL still
        # catches MONOTONE association, and a NON-MONOTONE hazard (events cluster at mid-
        # time: linear assoc ~0 but real dependence) is caught by the bias-corrected dCor
        # independence T-test (.dcor_t > 2, the 2-sigma universal). That test is n-INVARIANT
        # as a decision (a t-threshold, not total bits) and its U-centered statistic is ~0
        # under independence, so non-survival pairs do NOT flood as n grows -- and it only
        # fires AFTER the survival signature ss$pct>=30 gate above. Accept on EITHER route.
        if (sum(cc) <= 10) next
        tmc <- tm_vals[cc]; evc <- ev01[cc]
        assoc_bits <- .mdl_forward_gauss(tmc, matrix(evc, ncol = 1),
                                         fp = .col_fingerprint(evc))$bits_saved
        monotone <- is.finite(assoc_bits) && assoc_bits > 0
        nonmono  <- .dcor_t(evc, tmc) > 2                      # bias-corrected dCor independence reject
        if (!monotone && !nonmono) next                       # event independent of time -> not survival
        assoc <- abs(suppressWarnings(stats::cor(tmc, evc)))   # bounded ranking weight only
        if (is.na(assoc)) assoc <- 0
        evs <- c(evs, ev); tms <- c(tms, tm)
        ranks <- c(ranks, ss$score + 60 * assoc +
                     .name_bonus_pts(tm, name_bonus$time_variable) +
                     .name_bonus_pts(ev, name_bonus$event_variable))
        sss[[length(sss) + 1]] <- ss
      }
    }
  }
  if (length(evs) > 0) {
    # canonical pick: rank desc, then (event,time) value-fingerprint asc (C1).
    # NA-silence on an unbreakable tie (value-identical (event,time) pairs, e.g.
    # duplicate columns) -- no single pick is both relabel-invariant and
    # S_n-equivariant, so stay silent.
    fps <- unname(vapply(seq_along(evs), function(k)
      paste(.col_fingerprint(data[[evs[k]]]), .col_fingerprint(data[[tms[k]]]), sep = "~"),
      character(1L)))
    ord <- order(-ranks, fps)
    tie <- length(ord) >= 2L && isTRUE(ranks[ord[1]] == ranks[ord[2]]) &&
           isTRUE(fps[ord[1]] == fps[ord[2]])
    if (!tie) {
      w <- ord[1]
      best_event <- evs[w]; best_time <- tms[w]; best_ss <- sss[[w]]
      ss <- best_ss
      col <- if (component == "time") best_time else best_event
      return(list(
        found = TRUE, role = paste0(component, "_variable"), column = col,
        detected_by = "mathematical", score = ss$score, max_score = ss$max,
        pct = ss$pct, breakdown = ss$breakdown,
        event_col = best_event, time_col = best_time
      ))
    }
  }
  list(
    found = FALSE, role = paste0(component, "_variable"), column = NULL,
    detected_by = "mathematical", score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

# D3 ANCESTOR PORT (S_n-safe): long-format panel detection. A longitudinal panel has
# (1) a CLUSTER column with >= 3 distinct levels each REPEATED >= 3 times (so k <= n/3
# automatically -- no tuned upper cap); (2) a numeric TIME/visit column STRICTLY
# increasing WITHIN every cluster (the long-format ordering, and varying, not constant);
# (3) >= 1 numeric MEASURE column varying WITHIN clusters. This panel signal is what
# separates a subject-cluster (-> subject_id = cluster, repeated_measures = measures)
# from a cross-sectional grouping (a group_var has few levels but NO within-cluster
# monotone time and independent rows). All tests are column-VALUE functions (the cluster
# pick uses the value-fingerprint canonical key), so RELABEL/S_n hold; row order is a
# data property the column gates never permute. Returns NULL or list(cluster, time,
# measures). Strictly-increasing-within-cluster on clusters of size >= 3 makes a spurious
# panel on shuffled cross-sectional data vanishingly unlikely.
.detect_long_panel <- function(data, var_info) {
  n <- nrow(data); if (n < 20L) return(NULL)
  mono_inc_within <- function(cl, tv) {              # tv strictly increasing within EVERY cluster, and varies
    ok <- !is.na(cl) & is.finite(tv); cl <- cl[ok]; tv <- tv[ok]
    if (length(tv) < 3L) return(FALSE)
    varied <- FALSE
    for (lv in unique(cl)) {
      x <- tv[cl == lv]
      if (length(x) >= 3L) {
        if (any(diff(x) <= 0)) return(FALSE)        # not strictly increasing in stored order
        varied <- TRUE
      }
    }
    varied
  }
  within_var <- function(cl, mv) {                   # measure varies within clusters
    ok <- !is.na(cl) & is.finite(mv); cl <- cl[ok]; mv <- mv[ok]
    vs <- tapply(mv, cl, function(z) if (length(z) >= 2L) stats::var(z) else NA_real_)
    any(is.finite(vs) & vs > 0)
  }
  # cluster candidates: >= 3 levels, each repeated >= 3 times, not unique
  clust <- character(0)
  for (col in var_info$column) {
    v <- data[[col]]; nz <- v[!is.na(v)]
    if (length(nz) < 20L) next
    k <- length(unique(nz)); if (k < 3L || k >= length(nz)) next
    if (min(table(nz)) < 3L) next                    # every level repeated >= 3x => k <= n/3 (no tuned cap)
    clust <- c(clust, col)
  }
  if (length(clust) == 0L) return(NULL)
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  panels <- list()
  for (cc in clust) {
    cl <- as.character(data[[cc]])
    cand_num <- setdiff(num_cols, cc)
    tcols <- cand_num[vapply(cand_num, function(t) mono_inc_within(cl, suppressWarnings(as.numeric(data[[t]]))), logical(1L))]
    if (length(tcols) == 0L) next
    mcols <- cand_num[vapply(cand_num, function(m) within_var(cl, suppressWarnings(as.numeric(data[[m]]))), logical(1L))]
    mcols <- setdiff(mcols, tcols)
    if (length(mcols) < 1L) next
    panels[[length(panels) + 1L]] <- list(cluster = cc, time = tcols, measures = mcols)
  }
  if (length(panels) == 0L) return(NULL)
  # canonical pick among cluster candidates: smallest value-fingerprint (order-free)
  fps <- vapply(panels, function(p) .col_fingerprint(data[[p$cluster]]), character(1L))
  panels[[order(fps)[1]]]
}

.detect_subject_id <- function(data, var_info, name_bonus = NULL) {
  id_cols <- var_info$column[var_info$user_type == "ID"]
  id_pick <- if (length(id_cols) > 0) .canonical_pick(id_cols, rep(1, length(id_cols)), data) else NA_character_
  if (!is.na(id_pick)) {                                 # NA = unbreakable tie -> fall through
    return(list(
      found = TRUE, role = "subject_id", column = id_pick,
      detected_by = "mathematical", score = 90, max_score = 100, pct = 90,
      breakdown = list(list(name = "ID type", score = 90, max = 100, detail = "tagged as ID")),
      needs_confirmation = FALSE
    ))
  }
  id_like_cols <- character(0)
  for (col in var_info$column) {
    vals <- data[[col]]
    nonna <- vals[!is.na(vals)]                           # C2: NA-tolerant (was: any NA -> reject)
    if (length(nonna) < 10L) next
    if (length(unique(nonna)) != length(nonna)) next       # unique among OBSERVED values
    if (mean(is.na(vals)) > 0.2) next                      # too sparse to be a clean id
    id_like <- is.character(vals) || is.factor(vals) ||
      (is.numeric(vals) && isTRUE(all(nonna == round(nonna))))
    if (id_like) id_like_cols <- c(id_like_cols, col)
  }
  if (length(id_like_cols) > 0) {
    nb_hits <- id_like_cols[vapply(id_like_cols,
      function(col) .name_bonus_pts(col, name_bonus$subject_id) > 0, logical(1L))]
    pool <- if (length(nb_hits) > 0) nb_hits else id_like_cols
    pick <- .canonical_pick(pool, rep(1, length(pool)), data)           # order-free (C1)
    if (!is.na(pick)) return(list(
      found = TRUE, role = "subject_id", column = pick,
      detected_by = "mathematical", score = 70, max_score = 100, pct = 70,
      breakdown = list(list(
        name = "Unique ID-like values", score = 70, max = 100,
        detail = "all unique, integer/character (not numeric-continuous)"
      )),
      needs_confirmation = TRUE
    ))
  }
  # D3 cluster cross-check: in LONG-format data the subject id REPEATS (not unique), so
  # the unique-value paths above miss it. A subject CLUSTER is the cluster column of a
  # longitudinal panel (.detect_long_panel) -- 3-50 repeated levels WITH a within-cluster
  # monotone time + within-cluster-varying measures, which separates it from a group_var
  # (few levels, no panel structure). Value-only; name routes via name_bonus only.
  panel <- .detect_long_panel(data, var_info)
  if (!is.null(panel)) {
    return(list(
      found = TRUE, role = "subject_id", column = panel$cluster,
      detected_by = if (is.null(name_bonus$subject_id)) "mathematical" else "mathematical+name",
      score = 70, max_score = 100, pct = 70,
      breakdown = list(list(name = "Cluster id of a longitudinal panel", score = 70, max = 100,
                            detail = "3+ repeated levels with within-cluster monotone time + varying measures")),
      needs_confirmation = TRUE
    ))
  }
  list(
    found = FALSE, role = "subject_id", column = NULL, detected_by = "mathematical",
    score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

# v0.2.0 (phase G2): MDL collider/sink. An OUTCOME is a dependency SINK -- a column
# whose two-part code length drops below its marginal (null) code length when
# regressed on >= 2 predictors chosen by MDL FORWARD SELECTION (a predictor enters
# only if it pays for its (1/2)log2(n) parameter bits). This DISSOLVES the v0.1.x
# tuned cutoffs into one derived criterion:
#   CORR_MIN (which predictors)  -> a predictor enters iff it reduces total bits.
#   TAU      (how explained)     -> the n-aware MDL break-even (n/2)log2(1/(1-R^2)) > (k/2)log2(n).
#   GAMMA    (mutual independence)-> a collinear/redundant predictor adds NO new bits,
#       so a cluster member (repeated/scale/paired partner) selects ONE predictor and
#       is NOT flagged, while a true collider (>=2 independent causes) selects >=2.
# MIN_PRED=2 is the lone STRUCTURAL (definitional) constant; the bit margin is the
# only free calibration and is surfaced as confidence (G4). Causal DIRECTION stays
# unidentified under the linear-Gaussian symmetric sink (honest silence; G3 ANM
# breaks it where non-Gaussian/nonlinear). All features value-based => S_n/TURNUSOL
# hold; n<10 protection preserved (MDL never evades small-n -- (k/2)log2(n) is HARDER
# to pay at small n).
.detect_outcome <- function(data, var_info, type, name_bonus = NULL, claimed = character(0)) {
  # STRUCTURAL identifiability constants (NOT tuned): an OUTCOME is a dependency SINK.
  # CONTINUOUS needs >= 2 mutually-independent causes (y=a+b is symmetric: a=y-b looks
  # identical, so a lone cause cannot orient the sink). BINARY is asymmetric (a 0/1 column
  # is a function OF continuous features, not the reverse) so ONE cause already identifies
  # the effect -- BUT only once the columns a more-specific structural role already owns (a
  # survival EVENT, a SCALE-battery item, a repeated-measure, a pair partner) are removed
  # from contention. (DV regression triage: the v0.2.0 outcome detector lacked that
  # ATOMICITY and double-claimed those columns -- the survival event, whose one strong cause
  # is time, and collinear scale items out-ranked the true multi-cause outcome, dropping
  # outcome_binary below v0.1.x AND baseline. Root cause = missing atomicity, NOT MIN_PRED;
  # raising MIN_PRED to 2 for binary only crushed recall. Fix = `claimed` exclusion below.)
  MIN_PRED <- if (type == "continuous") 2L else 1L
  if (type == "continuous") {
    cands <- var_info$column[var_info$user_type == "Continuous"]
  } else {
    cands <- var_info$column[var_info$user_type %in% c("Binary", "Categorical")]
    cands <- cands[vapply(cands, function(c) length(unique(data[[c]][!is.na(data[[c]])])) == 2, logical(1L))]
  }
  # ATOMICITY (DV triage root-cause fix): an OUTCOME is not a column already owned by a
  # more-specific structural role (survival event/time, scale item, repeated measure, pair
  # partner, subject id, group). `claimed` is the union of those roles' columns detected
  # upstream (run_all_detections threads it in detection order). Value-identified columns =>
  # the exclusion is RELABEL / S_n invariant. This removes the event-steal and scale-steal
  # at the root, so binary keeps MIN_PRED=1 (recall) without the precision leak.
  cands <- setdiff(cands, claimed)
  none <- function(best = 0, reason = "no candidates")
    list(found = FALSE, role = paste0("outcome_", type), column = NULL,
         detected_by = "mathematical", score = 0, max_score = 100, pct = 0,
         breakdown = list(), confidence = round(best, 3), reason = reason)
  if (length(cands) == 0) return(none())
  numu <- var_info$column[var_info$user_type %in% c("Continuous", "Ordinal", "Binary")]
  if (length(numu) < 3) return(none(0, "too few columns for a collider"))

  bits  <- setNames(numeric(length(cands)), cands)   # MDL bits saved vs the null (marginal) model
  ksel  <- setNames(integer(length(cands)), cands)   # number of MDL-selected predictors
  indep <- setNames(rep(TRUE, length(cands)), cands)  # continuous: the selected causes must be mutually independent
  toppred <- setNames(rep(NA_character_, length(cands)), cands)  # top MDL predictor name (for ANM direction annotation)
  for (y in cands) {
    preds <- setdiff(numu, y)
    if (length(preds) < MIN_PRED) next
    yv <- suppressWarnings(as.numeric(data[[y]]))
    if (sum(is.finite(yv)) < 10L) next               # n<10 protection (inviolable)
    if (stats::sd(yv[is.finite(yv)]) < 1e-12) next
    X  <- vapply(preds, function(c) suppressWarnings(as.numeric(data[[c]])), numeric(nrow(data)))
    fp <- vapply(preds, function(c) .col_fingerprint(data[[c]]), character(1L))   # order-free tie-break
    # EVIDENCE MODEL by type (council D-C/D-D). CONTINUOUS: Gaussian MDL forward selection
    # (an LPM/RSS code is fine for a continuous sink). BINARY: a 0/1 sink is NON-Gaussian, so
    # the Gaussian-LPM RSS UNDERSTATES its logistic dependence (DV triage: outcome_binary fell
    # below v0.1.x AND the trivial type baseline). Score binary by the MI/dCor backbone (D-D)
    # instead -- count the causes the BIAS-CORRECTED dCor independence test rejects
    # (.dcor_t > 2, the 2-sigma universal). The U-centered .dcor_t (mean 0 under independence)
    # is used rather than raw .dep/.dcor, which carry a ~0.27 positive floor on a binary-vs-
    # continuous pair (dependency.R:81-89) that would flood FALSE causes through the council's
    # nominal floor. MIN_PRED=1 (binary asymmetric; atomicity already removed the event/scale
    # confounders) preserves recall; T>2 is n-INVARIANT so a weak small-n signal stays honestly
    # SILENT (Rissanen). No new tuned constant. (D-D's MIN_PRED=2 was overruled by its OWN
    # revert-bar: it collapses recall to ~0.28-0.33, below both baselines.)
    if (type == "continuous") {
      ccidx <- which(is.finite(yv) & (rowSums(!is.finite(X)) == 0))   # joint complete cases (design matrix)
      if (length(ccidx) < 10L) next
      if (length(ccidx) > .MDL_NMAX) ccidx <- ccidx[.dcor_idx(length(ccidx), .MDL_NMAX, seed = 7L)]
      yc <- yv[ccidx]; if (stats::sd(yc) < 1e-12) next
      Xc <- X[ccidx, , drop = FALSE]
      res <- .mdl_forward_gauss(yc, Xc, fp)
      bits[y] <- max(0, res$bits_saved); ksel[y] <- res$k
      if (res$k >= 1L) toppred[y] <- preds[res$sel[1]]          # dominant MDL cause (for ANM annotation)
      # collider identifiability guard: reject a CLUSTER source whose >=2 "causes" are
      # collinear (each MDL-explained by another) -- only an OUTCOME's causes are independent.
      if (res$k >= MIN_PRED && res$bits_saved > 0)
        indep[y] <- .mdl_set_independent(Xc[, res$sel, drop = FALSE], fp[res$sel])
    } else {
      # PAIRWISE association (NA-robust per predictor; .dcor_t row-caps + filters internally).
      # Joint complete-cases would drop the candidate whenever ANY predictor is character-
      # coerced-to-NA, starving recall -- so the binary path reads the FULL columns pairwise.
      assoc <- vapply(seq_along(preds), function(j) .dcor_t(yv, X[, j]), numeric(1L))
      paid  <- which(is.finite(assoc) & assoc > 2)              # 2-sigma bias-corrected dCor reject
      bits[y] <- if (length(paid)) sum(assoc[paid]) else 0       # aggregate association = value-based rank
      ksel[y] <- length(paid)
      if (length(paid)) toppred[y] <- preds[paid[which.max(assoc[paid])]]
    }
  }
  qual <- (ksel >= MIN_PRED) & (bits > 0) & indep    # collider: >=2 MDL-selected, net-compressing, mutually-independent causes
  best <- max(bits)
  if (!any(qual)) return(none(best, sprintf("no MDL collider sink (best bits saved=%.1f; need >=%d predictor(s))", best, MIN_PRED)))
  qcols <- cands[qual]
  kw <- name_bonus[[paste0("outcome_", type)]]
  pick <- .canonical_pick(qcols, unname(bits[qcols]), data)       # MATH pick: value-canonical; NA on identical tie
  if (is.na(pick)) return(none(best, "tied MDL collider candidates (silence)"))
  # Capped name hint (Tier-1, <= 10% weight): a name match overrides the math pick
  # ONLY when the named candidate's MDL bits are within 10% of the math-best, so the
  # name channel can never override a > 10% mathematical margin (the old pool-
  # restriction could). No-op under name_bonus = NULL (.name_bonus_pts returns 0).
  nb_hits <- qcols[vapply(qcols, function(c) .name_bonus_pts(c, kw) > 0, logical(1L))]
  if (length(nb_hits) > 0 && !(pick %in% nb_hits)) {
    cand <- .canonical_pick(nb_hits, unname(bits[nb_hits]), data)
    if (!is.na(cand) && bits[cand] >= bits[pick] * 0.9) pick <- cand
  }
  conf_bits <- unname(bits[pick]); n_pred <- unname(ksel[pick])
  # G4: confidence = the description-length MARGIN (bits) between the chosen outcome and
  # the next-best candidate (or the null model if it is the only one). The silence
  # threshold is DERIVED, not tuned: a margin below one parameter's code length
  # (1/2)log2(n) means the best and runner-up are within coding noise -> ambiguous ->
  # honest silence. (Replaces the bare exact-tie test with a bit-margin gate.)
  others <- setdiff(qcols, pick)
  runner <- if (length(others)) max(bits[others]) else 0
  margin <- conf_bits - runner
  n_eff  <- min(nrow(data), .MDL_NMAX)
  # Margin-silence applies to CONTINUOUS only: it resolves WHICH single continuous column
  # is the sink. Binary is exempt because a dataset legitimately has MULTIPLE binary sinks
  # that are DIFFERENT roles (the binary outcome AND the survival event are both predictable
  # 0/1 columns) -- silencing on their close bits would suppress a real outcome (it regressed
  # outcome_binary below the v0.1.1 dangerous floor). The margin is still surfaced as confidence.
  if (type == "continuous" && margin < .mdl_param_bits(1, n_eff))
    return(none(best, sprintf("ambiguous MDL margin (%.1f < %.1f bits) -> silence", margin, .mdl_param_bits(1, n_eff))))
  # D-B: ANM direction orientation, gated on NON-GAUSSIAN escape. ANM is identifiable ONLY
  # when the cause/effect pair is non-Gaussian or nonlinear (Hoyer 2009); in the linear-
  # Gaussian symmetric core (a=y-b looks identical) the direction is UNIDENTIFIABLE -> honest
  # SILENCE. So we attempt the orientation ONLY when `.gauss_margin` REJECTS joint normality
  # for the pair; when both columns are Gaussian-admissible we do NOT call ANM and the label
  # stays "collider-implied" (frozen behavior on the 69/69 linear-Gaussian datasets, bit-for-
  # bit). This is annotation-only (no column/role change), so it is F1-neutral and never
  # forces a pick. LiNGAM (non-Gaussian LINEAR) is DEFERRED: the dCor-residual statistic sits
  # below its noise floor (anm.R:22-24), so claiming it would fabricate unmeasurable coverage
  # (Gödel/Heisenberg/Hume) -- left as the documented SILENT limitation.
  direction <- "collider-implied (causes -> sink)"
  if (type == "continuous" && !is.na(toppred[pick]) &&
      !(.gauss_margin(data[[toppred[pick]]]) && .gauss_margin(data[[pick]]))) {
    ad <- .anm_direction(data[[toppred[pick]]], data[[pick]])
    direction <- if (identical(ad$direction, "x->y"))
                   sprintf("ANM %s -> %s (asym=%.3f)", toppred[pick], pick, ad$asym)
                 else "unidentifiable (residual-symmetric)"
  }
  pct <- round(100 * (1 - 2^(-margin / 8)))          # bounded display map of the bit-margin
  list(found = TRUE, role = paste0("outcome_", type), column = pick,
       detected_by = "mathematical", score = pct, max_score = 100, pct = pct,
       breakdown = list(list(name = if (type == "continuous") "MDL collider sink (>=2 independent causes)" else "MI/dCor binary sink (dCor-rejected causes)",
                             score = pct, max = 100,
                             detail = sprintf("evidence=%.1f, margin=%.1f, causes=%d; direction: %s", conf_bits, margin, n_pred, direction))),
       confidence = round(margin, 3), needs_confirmation = TRUE)
}

.detect_repeated_measures <- function(data, var_info) {
  # D3 LONG-format path FIRST (dedup): a longitudinal panel is a SPECIFIC signal
  # (cluster + within-cluster monotone time + within-cluster-varying measures) that the
  # WIDE correlation-block search below would mis-handle (it pulls the integer cluster id
  # into the block). When a panel is present the data is long-format, so the repeated
  # measurements are panel$measures (the within-subject-varying numeric columns, the time
  # axis and cluster id excluded). On cross-sectional (WIDE) data .detect_long_panel
  # returns NULL -> the WIDE search runs unchanged (a NO-OP on the seen corpus).
  panel <- .detect_long_panel(data, var_info)
  if (!is.null(panel)) {
    return(list(found = TRUE, role = "repeated_measures", column = panel$measures,
                detected_by = "mathematical", score = 70, max_score = 120, pct = 58.3,
                breakdown = list(list(name = "Long-format within-subject repeated measures",
                                      score = 70, max = 120,
                                      detail = sprintf("%d measure(s) varying within a longitudinal cluster panel",
                                                       length(panel$measures))))))
  }
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  if (length(num_cols) < 3) {
    return(list(found = FALSE, role = "repeated_measures", column = NULL,
                detected_by = "mathematical", score = 0, max_score = 120, pct = 0, breakdown = list()))
  }
  # v0.1.1 (audit Finding C1): order-free correlation-block search replaces the
  # contiguous-window scan (which could only find adjacent columns and was
  # position-dependent). Candidate blocks are built from the value-based
  # correlation structure, scored by the unchanged .score_repeated_signature
  # (pct >= 30 LEGIT-threshold preserved), and the winner is canonical.
  M <- as.matrix(data[, num_cols, drop = FALSE])
  # G1: block adjacency is distance correlation (nonlinear coupling), but a cheap
  # linear-correlation pre-filter (|r| >= 0.3) prunes the O(p^2) dCor adjacency to
  # a small candidate set per seed -- sound because repeated measurements are
  # linearly correlated by construction, so a linearly-independent column is not a
  # block member. The block is then DEFINED by dCor (>= 0.5). Both measures are
  # symmetric value functions => RELABEL / S_n invariant. Cuts wide-data cost ~20x.
  Clin <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  Clin[is.na(Clin)] <- 0
  # canonical (value-fingerprint) column order so .score_repeated_signature's
  # order-dependent "Trend" sub-score is computed on an order-INVARIANT sequence.
  fp_order <- function(cc) cc[order(vapply(cc, function(c) .col_fingerprint(data[[c]]), character(1L)))]
  blocks <- list()
  for (s in seq_along(num_cols)) {
    pre <- setdiff(which(abs(Clin[s, ]) >= 0.30), s)       # cheap linear pre-screen
    if (length(pre) == 0L) next
    dv <- vapply(pre, function(j) .dep(M[, s], M[, j]), numeric(1L))   # C1: completed signature .dep on survivors
    cand <- pre[dv >= 0.5]
    if (length(cand) == 0L) next
    cand <- cand[order(-dv[dv >= 0.5])]                     # strongest dCor first
    cand <- utils::head(unique(c(s, cand)), 12L)            # cap block size (v0.1.0 used ~10-window)
    members <- fp_order(num_cols[cand])                    # value-defined, order-free
    if (length(members) >= 3L) blocks[[length(blocks) + 1L]] <- members
  }
  if (length(blocks) > 0) blocks <- blocks[!duplicated(vapply(blocks, paste, character(1L), collapse = ""))]
  bcols <- list(); bscore <- numeric(0); brs <- list()
  for (b in blocks) {
    rs <- .score_repeated_signature(b, data)
    if (rs$pct >= 30) { bcols[[length(bcols) + 1L]] <- b; bscore <- c(bscore, rs$score); brs[[length(brs) + 1L]] <- rs }
  }
  if (length(bcols) > 0) {
    fps <- unname(vapply(bcols, function(b)
      paste(sort(vapply(b, function(c) .col_fingerprint(data[[c]]), character(1L))), collapse = "~"),
      character(1L)))
    ord <- order(-bscore, fps)
    tie <- length(ord) >= 2L && isTRUE(bscore[ord[1]] == bscore[ord[2]]) && isTRUE(fps[ord[1]] == fps[ord[2]])
    if (!tie) {                                            # NA-silence on value-identical block tie (C1)
      rs <- brs[[ord[1]]]
      return(list(found = TRUE, role = "repeated_measures", column = bcols[[ord[1]]],
                  detected_by = "mathematical", score = rs$score, max_score = rs$max,
                  pct = rs$pct, breakdown = rs$breakdown))
    }
  }
  list(found = FALSE, role = "repeated_measures", column = NULL,
       detected_by = "mathematical", score = 0, max_score = 120, pct = 0, breakdown = list())
}

.detect_scale_items <- function(data, var_info) {
  ord_cols <- var_info$column[var_info$user_type %in% c("Ordinal", "Continuous")]
  scale_candidates <- c()
  for (col in ord_cols) {
    vals <- data[[col]][!is.na(data[[col]])]
    if (length(vals) < 20) next
    if (!all(vals == round(vals))) next
    n_unique <- length(unique(vals))
    if (n_unique >= 3 && n_unique <= 10) scale_candidates <- c(scale_candidates, col)
  }
  if (length(scale_candidates) >= 5) {
    ss <- .score_scale_items_signature(scale_candidates, data)
    if (ss$pct >= 30) {
      return(list(found = TRUE, role = "scale_items", column = scale_candidates,
                  detected_by = "mathematical", score = ss$score, max_score = ss$max,
                  pct = ss$pct, breakdown = ss$breakdown, n_items = length(scale_candidates)))
    }
  }
  # D3 ANCESTOR PORT (complementary): a DICHOTOMOUS IRT battery -- >= 8 binary {0,1}
  # columns (value-only; the ordinal path above misses binaries). Only when the ordinal
  # path yielded no scale, so it never overrides the strong ordinal detector. STRUCTURAL
  # precision guard with NO tuned threshold: reject a one-hot PARTITION (every respondent
  # endorses exactly one item -> max row-sum <= 1), which is an encoding, not a scale; a
  # real battery has respondents endorsing several items (max row-sum > 1).
  if (length(scale_candidates) < 5) {
    bin_cols <- character(0)
    for (col in var_info$column) {
      v <- data[[col]]; nz <- v[!is.na(v)]
      if (length(nz) < 20L || !is.numeric(v)) next
      u <- sort(unique(nz), method = "radix")
      if (length(u) == 2L && isTRUE(all(abs(u - c(0, 1)) < 1e-9))) bin_cols <- c(bin_cols, col)
    }
    if (length(bin_cols) >= 8L) {
      mb <- as.matrix(data[, bin_cols, drop = FALSE]); mb <- mb[stats::complete.cases(mb), , drop = FALSE]
      if (nrow(mb) >= 20L && max(rowSums(mb)) > 1) {
        ss <- .score_scale_items_signature(bin_cols, data)
        if (ss$pct >= 30) {
          return(list(found = TRUE, role = "scale_items", column = bin_cols,
                      detected_by = "mathematical", score = ss$score, max_score = ss$max,
                      pct = ss$pct, breakdown = ss$breakdown, n_items = length(bin_cols)))
        }
      }
    }
  }
  list(found = FALSE, role = "scale_items", column = NULL,
       detected_by = "mathematical", score = 0, max_score = 100, pct = 0, breakdown = list())
}

.detect_covariates <- function(data, var_info) {
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  cat_cols <- var_info$column[var_info$user_type %in% c("Categorical", "Binary")]
  covs <- c(num_cols, cat_cols)
  if (length(covs) > 0) {
    return(list(found = TRUE, role = "covariate", column = covs,
                detected_by = "mathematical", score = 50, max_score = 100, pct = 50,
                breakdown = list()))
  }
  list(found = FALSE, role = "covariate", column = NULL,
       detected_by = "mathematical", score = 0, max_score = 100, pct = 0, breakdown = list())
}

# Flatten a role-detection result to its claimed column names (scalar, vector, or the
# list(list(col1, col2)) pair shape). Value-identified, so the atomicity exclusion built
# from these is RELABEL / S_n invariant.
.role_cols <- function(r) {
  col <- r$column
  if (is.null(col)) return(character(0))
  if (is.list(col)) return(unname(unlist(lapply(col, function(p) unlist(p)))))
  as.character(col)
}

# Dispatch one role need to its detector. `claimed` = columns already owned by a more-
# specific structural role (threaded by run_all_detections); only the outcome detectors
# consult it (atomicity).
detect_variable_role <- function(data, var_info, need, name_bonus = NULL, claimed = character(0)) {
  switch(need,
    group_var          = .detect_group_var(data, var_info, name_bonus),
    paired_pairs       = .detect_paired_pairs(data, var_info),
    agreement_pairs    = .detect_agreement_pairs(data, var_info),
    time_variable      = .detect_survival_components(data, var_info, "time", name_bonus),
    event_variable     = .detect_survival_components(data, var_info, "event", name_bonus),
    subject_id         = .detect_subject_id(data, var_info, name_bonus),
    outcome_continuous = .detect_outcome(data, var_info, "continuous", name_bonus, claimed),
    outcome_binary     = .detect_outcome(data, var_info, "binary", name_bonus, claimed),
    repeated_measures  = .detect_repeated_measures(data, var_info),
    scale_items        = .detect_scale_items(data, var_info),
    covariate          = .detect_covariates(data, var_info),
    list(found = FALSE, role = need, column = NULL, detected_by = "none",
         score = 0, max_score = 0, pct = 0, breakdown = list())
  )
}

# The canonical role-need set, in detection order. Structural roles (survival, subject id,
# repeated, scale) resolve BEFORE the outcomes so their columns are `claimed` (atomicity).
.ROLE_NEEDS <- c(
  "group_var", "paired_pairs", "agreement_pairs", "time_variable",
  "event_variable", "subject_id", "repeated_measures", "scale_items",
  "outcome_continuous", "outcome_binary", "covariate"
)

# Roles whose detected columns BLOCK a column from also being an outcome (atomicity).
# Excludes covariate (the catch-all) and the outcome roles themselves.
.CLAIMING_ROLES <- c("group_var", "paired_pairs", "agreement_pairs", "time_variable",
                     "event_variable", "subject_id", "repeated_measures", "scale_items")

# Run every detector + per-column value classification + potential pairs.
run_all_detections <- function(data, var_info, name_bonus = NULL, verbose = FALSE) {
  roles <- list()
  for (need in .ROLE_NEEDS) {
    claimed <- unique(unlist(lapply(intersect(names(roles), .CLAIMING_ROLES),
                                    function(rn) .role_cols(roles[[rn]]))))
    roles[[need]] <- tryCatch(
      detect_variable_role(data, var_info, need, name_bonus, claimed),
      error = function(e) list(found = FALSE, role = need, .error = conditionMessage(e))
    )
    .say(verbose, "[rolescry] %s: found=%s", need, isTRUE(roles[[need]]$found))
  }
  value_types <- vapply(names(data), function(col) {
    tryCatch(classify_value_type(data[[col]])$type, error = function(e) NA_character_)
  }, character(1L))
  potential_pairs <- list()
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  if (length(num_cols) >= 2) {
    for (i in seq_along(num_cols)[-length(num_cols)]) {
      for (j in (i + 1):length(num_cols)) {
        a <- data[[num_cols[i]]]
        b <- data[[num_cols[j]]]
        ps <- .score_paired_signature(a, b)
        as <- .score_agreement_signature(a, b)
        if (ps$pct >= 30 || as$pct >= 30) {
          potential_pairs[[length(potential_pairs) + 1]] <- list(
            col1 = num_cols[i], col2 = num_cols[j],
            paired_score = ps, agreement_score = as
          )
        }
      }
    }
  }
  list(roles = roles, value_types = value_types, potential_pairs = potential_pairs)
}
