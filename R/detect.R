# Role detection engine. Ported from MDStatR engine/run_097_detect_roles.R
# (detect_variable_role + .detect_* + run_all_detections), with the Shiny rv
# coupling severed: var_info is an explicit data.frame argument, and the
# confirmed_roles / config layers (Shiny-runtime concerns) are dropped -- the
# Phase A baseline ran with config = NULL and confirmed_roles = NULL, so it
# already entered the mathematical layer directly, making this exactly
# equivalent. The optional capped name bonus is threaded into the
# candidate-selecting detectors (group, outcome, subject id, survival); with
# name_bonus = NULL it is a no-op (purely mathematical). No exports.
#
# State at v0.3.0 (pass 5) -- two layers, read apart:
#   * DERIVED deciders (no point scores in the decision): .detect_group_var (normalized
#     entropy), .detect_pairs_generic (paired vs agreement by the MDL shift sign
#     .mdl_shift_bits, pair evidence .mdl_pair_bits + dCor fallback), ordinal scale_items
#     (.mdl_factor_bits one-factor code), .detect_outcome (.mdl_forward_gauss / Bernoulli
#     MDL / .dcor_t, .anm_direction orientation), the long-format panel path, coupled_pair.
#   * PORTED point scorers that REMAIN (signatures.R, v0.1.0 logic, pct >= 30 accept):
#     .detect_survival_components (signature gate), wide .detect_repeated_measures, and the
#     binary (IRT) scale path. They remain because their derived replacements were measured
#     and rejected under the no-regression rule (see the v0.3.0 development report, rolescry_v030/v030_dev_report.md, not distributed with the package).
# Origin: see "Relationship to MDStatR" in README.md

# Normalized Shannon entropy (bits) of a categorical column: B = H / log2(k), B = 1 iff
# perfectly balanced. Value-only; NA dropped. NULL when k < 2 or any level seen < 3 times.
.group_balance <- function(x) {
  cnt <- as.integer(table(x[!is.na(x)]))
  cnt <- cnt[cnt > 0L]
  k <- length(cnt)
  if (k < 2L || any(cnt < 3L)) return(NULL)
  p <- cnt / sum(cnt)
  list(B = round(min(1, -sum(p * log2(p)) / log2(k)), 12), k = k)   # round: float-noise guard for exact rank ties, not a tuned constant
}

.detect_group_var <- function(data, var_info, name_bonus = NULL) {
  # Rule (pass 5): group_var = the char/factor Categorical/Binary column with the highest
  # normalized entropy B (+ capped name bonus); ties -> fewer levels -> value fingerprint;
  # value-identical top two -> silence (as .canonical_pick). Candidates stay char/factor-only:
  # a numeric-coded treatment is indistinguishable name-blind from a low-card covariate/
  # outcome (Lemma-2 wall). Measured-and-rejected alternatives: see the v0.3.0 development report (rolescry_v030/v030_dev_report.md, not distributed with the package).
  candidates <- var_info$column[var_info$user_type %in% c("Categorical", "Binary")]
  cols <- character(0); ranks <- numeric(0); gbs <- list()
  for (col in candidates) {
    x <- data[[col]]
    if (!is.character(x) && !is.factor(x)) next
    gb <- .group_balance(x)
    if (is.null(gb)) next
    cols <- c(cols, col)
    ranks <- c(ranks, 100 * gb$B + .name_bonus_pts(col, name_bonus$group_var))
    gbs[[length(gbs) + 1]] <- gb
  }
  pick <- NA_character_
  if (length(cols) > 0) {
    ks  <- vapply(gbs, function(g) g$k, integer(1L))
    fps <- unname(vapply(cols, function(k) .col_fingerprint(data[[k]]), character(1L)))
    ord <- order(-ranks, ks, fps)
    i1 <- ord[1]; i2 <- ord[2]
    tie <- length(ord) >= 2L && ranks[i1] == ranks[i2] && ks[i1] == ks[i2] && fps[i1] == fps[i2]
    if (!tie) pick <- cols[i1]                           # NA = unbreakable tie -> silence (C1)
  }
  if (!is.na(pick)) {
    gb <- gbs[[match(pick, cols)]]
    pct <- round(100 * gb$B, 1)
    return(list(
      found = TRUE, role = "group_var", column = pick,
      detected_by = if (is.null(name_bonus$group_var)) "mathematical" else "mathematical+name",
      score = pct, max_score = 100, pct = pct,
      breakdown = list(list(name = "normalized entropy", score = gb$B, max = 1,
                            detail = sprintf("%d levels, H/log2k=%.3f", gb$k, gb$B)))
    ))
  }
  list(
    found = FALSE, role = "group_var", column = NULL, detected_by = "mathematical",
    score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

# Shared admission rule of the pair detectors (also used by verify_roles): the
# consensus+difference MDL closed form, with the dCor fallback for non-Gaussian margins.
.pair_admitted <- function(a, b) {
  mb <- .mdl_pair_bits(a, b)
  cf <- isTRUE(mb$bits > 0 && mb$agree)
  accept <- cf
  route <- if (cf) "closed-form" else "none"
  if (!cf && !(.gauss_margin(a) && .gauss_margin(b))) {
    m_eff <- sum(is.finite(a) & is.finite(b))
    accept <- isTRUE(mb$agree) && .dep(a, b) > 2 / sqrt(max(m_eff, 1))
    if (accept) route <- "dcor-fallback"
  }
  list(admitted = accept, bits = mb$bits, agree = mb$agree, route = route)
}

.detect_pairs_generic <- function(data, var_info, role) {
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  if (length(num_cols) < 2) {
    return(list(found = FALSE, role = role, column = NULL, detected_by = "mathematical",
                score = 0, max_score = NA, pct = 0, breakdown = list()))
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
      # .mdl_pair_bits test is O(n), so it runs on all pairs cheaply; no point scorer
      # is involved in the pair decision.
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
      pa <- .pair_admitted(a, b)
      mb <- list(bits = pa$bits)
      accept <- pa$admitted
      if (accept) {
        pairs[[length(pairs) + 1]] <- list(col1 = num_cols[i], col2 = num_cols[j],
                                           bits = mb$bits)
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
      bits <- top$bits
      # the single best pair is classified ONCE: systematic shift => paired, else agreement
      shift <- .mdl_shift_bits(data[[top$col1]], data[[top$col2]])
      if ((role == "paired_pairs") != (shift > 0)) {
        return(list(found = FALSE, role = role, column = NULL, detected_by = "mathematical",
                    score = 0, max_score = NA, pct = 0, breakdown = list(),
                    reason = sprintf("pair reported under %s",
                                     if (role == "paired_pairs") "agreement_pairs" else "paired_pairs")))
      }
      return(list(
        found = TRUE, role = role,
        column = list(list(col1 = top$col1, col2 = top$col2)),
        all_pairs = pairs[ord], detected_by = "mathematical",
        score = bits, max_score = NA, pct = max(0, round(100 * (1 - 2^(-bits / 8)))),  # display map, not a decision
        breakdown = list(list(name = "consensus/difference MDL", score = bits, max = NA,
                              detail = sprintf("pair bits=%.1f, shift bits=%.1f", bits, shift)))
      ))
    }
  }
  list(found = FALSE, role = role, column = NULL, detected_by = "mathematical",
       score = 0, max_score = NA, pct = 0, breakdown = list())
}

.detect_paired_pairs <- function(data, var_info) {
  .detect_pairs_generic(data, var_info, "paired_pairs")
}
.detect_agreement_pairs <- function(data, var_info) {
  .detect_pairs_generic(data, var_info, "agreement_pairs")
}

# 0.3.0: conditional identification of the SURVIVAL reading of a dependent
# (binary, positive-continuous) pair. The identifiable object from P alone is the
# PAIR; whether it is (time, event) or (cause, binary outcome) is design-relative
# (same ceiling class as Lemma 2). Censoring geometry singles out the survival
# reading when ALL of: (i) the continuous member is anchored at zero (mean - 2*sd
# <= 0, the 2-sigma universal -- a duration is ratio-scale with 0 in reach, an
# interval-scale measurement is not); (ii) events sit strictly earlier than
# censorings (median); (iii) the top time-tertile's event rate is strictly below
# the bottom's -- under censoring the longest observed times are mostly censored.
# All tests are value-only => RELABEL/S_n invariant. Documented ceiling: a
# zero-anchored continuous CAUSE whose larger values are protective mimics this
# geometry and still reads as survival.
.censoring_decisive <- function(tv, ev) {
  cc <- is.finite(tv) & !is.na(ev)
  tv <- tv[cc]; ev <- ev[cc]
  if (length(tv) < 20L || length(unique(ev)) != 2L) return(FALSE)   # n<20 estimator floor
  if (!any(ev == 0) || !any(ev == 1)) return(FALSE)                 # no censoring observed -> geometry undefined
  sdv <- stats::sd(tv)
  if (!is.finite(sdv) || sdv <= 0) return(FALSE)
  # (i) zero-anchored scale: a duration is a ratio-scale quantity anchored at 0,
  # so 0 lies inside its 2-sigma envelope (mean - 2*sd <= 0). An interval-scale
  # measurement far from zero (blood pressure, albumin) fails this. Reuses the
  # package's 2-sigma universal; no new constant.
  if (mean(tv) - 2 * sdv > 0) return(FALSE)
  # (ii) events sit strictly earlier than censorings
  if (!(stats::median(tv[ev == 1]) < stats::median(tv[ev == 0]))) return(FALSE)
  # (iii) the top time-tertile's event rate is strictly below the bottom's:
  # under censoring, the longest observed times are mostly censored. Only the
  # ENDS are compared -- the middle tertile follows the hazard shape and is not
  # constrained by censoring.
  q <- stats::quantile(tv, c(1/3, 2/3), names = FALSE, type = 1)    # order statistics: value-only
  g <- findInterval(tv, q)
  r <- vapply(c(0, 2), function(k) if (any(g == k)) mean(ev[g == k]) else NA_real_, numeric(1))
  if (anyNA(r)) return(FALSE)
  r[1] > r[2]
}

# Survival association gate shared by the survival detector and verify_roles: event
# compresses time under a 1-predictor Gaussian MDL model (monotone) OR the bias-corrected
# dCor independence T-test rejects (non-monotone); evaluated on complete cases.
.survival_assoc <- function(tv, ev01, min_n = 10L) {
  cc <- !is.na(tv) & !is.na(ev01)
  if (sum(cc) <= min_n) return(list(admitted = FALSE, bits = NA_real_, t = 0, monotone = FALSE, nonmono = FALSE))
  tvc <- tv[cc]; evc <- ev01[cc]
  bits <- .mdl_forward_gauss(tvc, matrix(evc, ncol = 1), fp = .col_fingerprint(evc))$bits_saved
  monotone <- is.finite(bits) && bits > 0
  tt <- .dcor_t(evc, tvc)
  nonmono <- tt > 2
  list(admitted = monotone || nonmono, bits = bits, t = tt, monotone = monotone, nonmono = nonmono)
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
  # 0.3.0: an all-distinct positive integer column value-types as ID, but a follow-up
  # time recorded in whole days has exactly that signature (ovarian futime). Admit
  # numeric ID-typed columns as TIME candidates; the event-association gate below
  # decides, so a true identifier (independent of every event) is not admitted.
  id_time <- var_info$column[var_info$user_type == "ID" &
                vapply(var_info$column, function(c) is.numeric(data[[c]]), logical(1L))]
  time_cands <- c(numeric_cols, id_time)
  evs <- character(0); tms <- character(0); ranks <- numeric(0); sss <- list(); nbs <- numeric(0)
  acs <- numeric(0); nns <- integer(0)                 # |cor| + n per admitted pair (dominance gate)
  for (ev in event_cands) {
    ev01 <- .normalize_event(data[[ev]])                 # C3: {1,2}/factor/logical/{0,1,2} -> {0,1}
    if (is.null(ev01)) next
    obs <- !is.na(ev01)                                  # C2: available-case, no hard NA skip
    if (sum(obs) < 10) next                              # n<10 LEGIT-threshold preserved
    if (length(unique(ev01[obs])) != 2L) next
    for (tm in time_cands) {
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
        # Survival association, two routes: the linear-MDL break-even floor (event compresses
        # time under a 1-predictor Gaussian model) catches MONOTONE association; the
        # bias-corrected dCor independence T-test (.dcor_t > 2, the 2-sigma universal; n-invariant,
        # ~0 under independence) catches NON-MONOTONE hazard. Accept on EITHER route, only after the
        # survival signature ss$pct>=30 gate above. Measured-and-rejected alternatives: see the v0.3.0 development report (rolescry_v030/v030_dev_report.md, not distributed with the package).
        if (sum(cc) <= 10) next
        tmc <- tm_vals[cc]; evc <- ev01[cc]
        sa <- .survival_assoc(tm_vals, ev01)
        monotone <- sa$monotone; nonmono <- sa$nonmono
        if (!sa$admitted) next                                # event independent of time -> not survival
        # An ID-typed time candidate overturns a structural typing, so it must pass BOTH
        # evidence routes (monotone MDL AND dCor independence reject), not either one;
        # a random identifier passes one route by chance far more often than both.
        if (tm %in% id_time && !(monotone && nonmono)) next
        assoc <- abs(suppressWarnings(stats::cor(tmc, evc)))   # bounded ranking weight only
        if (is.na(assoc)) assoc <- 0
        evs <- c(evs, ev); tms <- c(tms, tm)
        acs <- c(acs, assoc); nns <- c(nns, sum(cc))
        nbs <- c(nbs, .name_bonus_pts(tm, name_bonus$time_variable) +
                        .name_bonus_pts(ev, name_bonus$event_variable))
        ranks <- c(ranks, ss$score + 60 * assoc + nbs[length(nbs)])
        sss[[length(sss) + 1]] <- ss
      }
    }
  }
  if (length(evs) > 1L) {
    # 0.3.0: re-evaluate every admitted pair on ONE common row set (rows complete in all
    # admitted candidate columns). Each pair is first judged on its own available cases,
    # so pairs whose columns are missing in different rows were compared on different
    # subsamples and one could win on noise (lung: sex/meal.cal over status/time).
    # Gates and ranking are unchanged; only the rows they are computed on are common.
    common <- stats::complete.cases(data[, unique(c(evs, tms)), drop = FALSE])
    if (sum(common) >= 10L && sum(!common) > 0L) {
      dc <- data[common, , drop = FALSE]; keep <- logical(length(evs))
      for (k in seq_along(evs)) {
        e <- .normalize_event(dc[[evs[k]]]); t <- suppressWarnings(as.numeric(dc[[tms[k]]]))
        ss <- .score_survival_signature(tms[k], evs[k], dc)
        if (is.null(e) || length(unique(e)) != 2L || ss$pct < 30) next
        sa <- .survival_assoc(t, e, min_n = 9L)            # common block has no n<=10 skip (>= 10 rows)
        mono <- sa$monotone; nonm <- sa$nonmono
        if (!sa$admitted) next
        if (tms[k] %in% id_time && !(mono && nonm)) next
        a <- abs(suppressWarnings(stats::cor(t, e))); if (is.na(a)) a <- 0
        ranks[k] <- ss$score + 60 * a + nbs[k]; sss[[k]] <- ss; keep[k] <- TRUE
        acs[k] <- a; nns[k] <- sum(is.finite(t) & !is.na(e))
      }
      evs <- evs[keep]; tms <- tms[keep]; ranks <- ranks[keep]; sss <- sss[keep]; nbs <- nbs[keep]
      acs <- acs[keep]; nns <- nns[keep]
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
      # 0.3.0: assign the survival roles only when censoring geometry singles the
      # survival reading out; otherwise the honest output is the PAIR itself
      # (component "pair" -> the coupled_pair role), with the interpretation left
      # to the user. See .censoring_decisive.
      decisive <- .censoring_decisive(suppressWarnings(as.numeric(data[[best_time]])),
                                      .normalize_event(data[[best_event]]))
      # 0.3.0 dominance gate: Lemma 2 makes a SINGLE decisive pair's geometry
      # non-discriminating between readings, but when SEVERAL pairs pass the
      # decisive test the data offers competing survival readings. The winner is
      # assigned the survival roles only when its event-time association DOMINATES
      # every decisive competitor at the 2-sigma universal (Fisher z two-correlation
      # comparison, the same device as the panel-axis gate); otherwise the honest
      # output is the pair itself with the competition named.
      rival <- NA_character_
      if (decisive && length(evs) > 1L) {
        others <- setdiff(seq_along(evs), w)
        # Only DISJOINT-time rivals compete as READINGS: a rival sharing the winner's
        # time column agrees the table is survival-shaped on that axis and merely
        # offers another event candidate -- a ranking question the association weight
        # and the common-row re-evaluation already answer (pbc: time/hepato must not
        # veto time/status). A rival on a DIFFERENT time axis disputes which column
        # is the duration at all -- that is the real underdetermination.
        others <- others[tms[others] != tms[w]]
        # An ID-typed time candidate enters on an exceptional evidence bar (both
        # routes) because its prior reading is "identifier"; it may WIN on that
        # evidence but may not VETO a structurally-typed reading (pbc: a sequential
        # registry id correlates with administrative censoring and mimics survival).
        if (!(tms[w] %in% id_time)) others <- others[!(tms[others] %in% id_time)]
        odec <- others[vapply(others, function(k) .censoring_decisive(
          suppressWarnings(as.numeric(data[[tms[k]]])), .normalize_event(data[[evs[k]]])), logical(1L))]
        if (length(odec)) {
          cpt <- odec[which.max(ranks[odec])]
          zw <- atanh(min(acs[w], 1 - 1e-12)); zc <- atanh(min(acs[cpt], 1 - 1e-12))
          dominant <- nns[w] > 3L && nns[cpt] > 3L &&
            (zw - zc) / sqrt(1 / (nns[w] - 3) + 1 / (nns[cpt] - 3)) > 2
          if (!dominant) { decisive <- FALSE; rival <- paste(tms[cpt], evs[cpt], sep = "/") }
        }
      }
      if (component == "pair") {
        if (!decisive) return(list(
          found = TRUE, role = "coupled_pair",
          column = c(best_time, best_event),
          continuous_col = best_time, binary_col = best_event,
          interpretations = c(
            "survival: time = continuous member, event = binary member",
            "causal: cause = continuous member, binary outcome = binary member"),
          detected_by = "mathematical", score = ss$score, max_score = ss$max,
          pct = ss$pct, breakdown = ss$breakdown,
          reason = if (is.na(rival))
            "dependent binary/positive-continuous pair; censoring geometry not decisive, so the survival-vs-outcome reading is design-relative"
          else sprintf("several admissible survival readings without a 2-sigma dominance margin (competing pair: %s); the reading is left to the user", rival)
        ))
      } else if (decisive) {
        col <- if (component == "time") best_time else best_event
        return(list(
          found = TRUE, role = paste0(component, "_variable"), column = col,
          detected_by = "mathematical", score = ss$score, max_score = ss$max,
          pct = ss$pct, breakdown = ss$breakdown,
          event_col = best_event, time_col = best_time
        ))
      }
    }
  }
  list(
    found = FALSE,
    role = if (component == "pair") "coupled_pair" else paste0(component, "_variable"),
    column = NULL,
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
  dt_cols  <- var_info$column[var_info$user_type == "Datetime"]   # 0.3.0: a calendar column is a panel time axis
  panels <- list()
  for (cc in clust) {
    cl <- as.character(data[[cc]])
    cand_num <- setdiff(num_cols, cc)
    cand_t <- c(cand_num, setdiff(dt_cols, cc))
    tcols <- cand_t[vapply(cand_t, function(t) mono_inc_within(cl, .datetime_num(data[[t]])), logical(1L))]
    if (length(tcols) == 0L) next
    mcols <- cand_num[vapply(cand_num, function(m) within_var(cl, suppressWarnings(as.numeric(data[[m]]))), logical(1L))]
    mcols <- setdiff(mcols, tcols)
    if (length(mcols) < 1L) next
    panels[[length(panels) + 1L]] <- list(cluster = cc, time = tcols, measures = mcols)
  }
  if (length(panels) == 0L) return(NULL)
  if (length(panels) == 1L) return(panels[[1L]])
  # 0.3.0: a BALANCED panel is value-symmetric -- (firm, month) and (month, firm)
  # both satisfy the structural test, so the entity/time reading is not a function
  # of P alone. Discriminator: along the TRUE time axis the measures are serially
  # dependent within the entity. Pooled within-cluster lag-1 correlation, tested at
  # the 2-sigma universal (|r| sqrt(n) / sqrt(1 - r^2) > 2). Supported readings win;
  # several or none supported and the clusters differ -> AMBIGUOUS, reported as such
  # (the coupled_pair philosophy: name the symmetry, let the user decide).
  lag1_stats <- function(p) {
    cl <- as.character(data[[p$cluster]]); tv <- .datetime_num(data[[p$time[1]]])
    best <- list(t = 0, z = 0, w = 0)
    for (m in p$measures) {
      mv <- suppressWarnings(as.numeric(data[[m]]))
      # Two-way within transformation (panel fixed effects): subtract cluster and
      # time-position means. Raw pooling lets between-cluster level variance and
      # FIXED cross-sectional patterns (the same adjacent-position pairing recurring
      # in every cluster) masquerade as serial dependence on BOTH axis readings;
      # after demeaning only genuine serial structure survives.
      ok <- is.finite(mv) & is.finite(tv) & !is.na(cl)
      mv[ok] <- mv[ok] - stats::ave(mv[ok], cl[ok]) - stats::ave(mv[ok], tv[ok]) + mean(mv[ok])
      a <- c(); b <- c()
      for (lv in unique(cl)) {
        i <- which(cl == lv & is.finite(tv) & is.finite(mv))
        if (length(i) < 3L) next
        i <- i[order(tv[i])]
        a <- c(a, mv[i][-length(i)]); b <- c(b, mv[i][-1])
      }
      if (length(a) >= 8L && stats::sd(a) > 0 && stats::sd(b) > 0) {
        r <- stats::cor(a, b)
        if (is.finite(r) && abs(r) < 1) {
          tt <- abs(r) * sqrt(length(a)) / sqrt(1 - r^2)
          if (tt > best$t) best <- list(t = tt, z = atanh(abs(r)), w = length(a) - 3L)
        }
      }
    }
    best
  }
  st <- lapply(panels, lag1_stats)
  ts <- vapply(st, `[[`, numeric(1L), "t")
  sup <- which(ts > 2)
  if (length(sup) > 1L) {
    # Both readings pass the marginal gate. On the FALSE axis the recurring
    # adjacent-position pairing makes products dependent across clusters, so its
    # nominal t overstates the evidence; the TRUE axis carries structural serial
    # correlation bounded away from zero. Fisher z two-correlation comparison at
    # the same 2-sigma universal: a decisive gap singles the reading out,
    # otherwise the symmetry is real and stays reported as ambiguous.
    o <- sup[order(-ts[sup])]
    z1 <- st[[o[1]]]; z2 <- st[[o[2]]]
    if (z1$w > 0 && z2$w > 0 &&
        (z1$z - z2$z) / sqrt(1 / z1$w + 1 / z2$w) > 2) sup <- o[1]
  }
  pool <- if (length(sup) >= 1L) panels[sup] else panels
  fps <- vapply(pool, function(p) .col_fingerprint(data[[p$cluster]]), character(1L))
  if (length(unique(vapply(pool, `[[`, "", "cluster"))) == 1L || length(sup) == 1L)
    return(pool[[order(fps)[1]]])
  list(ambiguous = TRUE, panels = pool[order(fps)])
}

# Identifier-like column: >= 10 observed, unique among observed, <= 20 % missing, and
# character/factor or integer-valued numeric. Shared by the detector and verify_roles.
.id_like <- function(x) {
  nonna <- x[!is.na(x)]                                   # C2: NA-tolerant (was: any NA -> reject)
  if (length(nonna) < 10L) return(FALSE)
  if (length(unique(nonna)) != length(nonna)) return(FALSE)
  if (mean(is.na(x)) > 0.2) return(FALSE)
  is.character(x) || is.factor(x) || (is.numeric(x) && isTRUE(all(nonna == round(nonna))))
}

.detect_subject_id <- function(data, var_info, name_bonus = NULL, claimed = character(0)) {
  # 0.3.0: a column already owned by a survival role (e.g. an all-distinct integer
  # follow-up time) is not also a subject identifier (atomicity).
  id_cols <- setdiff(var_info$column[var_info$user_type == "ID"], claimed)
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
  for (col in setdiff(var_info$column, claimed)) {
    if (.id_like(data[[col]])) id_like_cols <- c(id_like_cols, col)
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
  if (!is.null(panel) && isTRUE(panel$ambiguous)) {
    # 0.3.0: axis-ambiguous balanced panel -> no single entity column is
    # identifiable; stay silent (the repeated_measures entry carries the readings).
    return(list(found = FALSE, role = "subject_id", column = NULL,
                detected_by = "mathematical", score = 0, max_score = 100, pct = 0,
                breakdown = list(),
                reason = "balanced panel: entity/time axis value-symmetric; see repeated_measures$panel_axes"))
  }
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

  bits  <- stats::setNames(numeric(length(cands)), cands)   # MDL bits saved vs the null (marginal) model
  ksel  <- stats::setNames(integer(length(cands)), cands)   # number of MDL-selected predictors
  indep <- stats::setNames(rep(TRUE, length(cands)), cands)  # continuous: the selected causes must be mutually independent
  toppred <- stats::setNames(rep(NA_character_, length(cands)), cands)  # top MDL predictor name (for ANM direction annotation)
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
    # SILENT (Rissanen). No new tuned constant. Measured-and-rejected alternatives: see the v0.3.0 development
    # report (rolescry_v030/v030_dev_report.md, not distributed with the package).
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
      # Bernoulli marginal MDL per predictor (NA-robust pairwise, as before: joint
      # complete-cases would starve recall when predictors are missing in different
      # rows). Replaces the bare 2-sigma dCor test, whose m-way scan inflated the
      # phantom-outcome rate on pure noise to ~1-(1-0.023)^m (~14% at m=4-5) and
      # whose familywise correction cost more true detections than it removed false
      # ones. The MDL code is selection-consistent by construction: each admitted
      # cause must pay its slope's 1/2 log2(n) parameter bits PLUS the log2(m)
      # model-index bits naming it among the m candidates (.mdl_bern_bits).
      yb <- data[[y]]; lv <- sort(unique(yb[!is.na(yb)]))
      y01 <- ifelse(is.na(yb), NA_integer_, as.integer(yb == lv[2L]))
      bb <- vapply(seq_along(preds), function(j)
        .mdl_bern_bits(y01, X[, j], length(preds)), numeric(1L))
      paid  <- which(is.finite(bb) & bb > 0)                    # net-compressing causes
      bits[y] <- if (length(paid)) sum(bb[paid]) else 0
      ksel[y] <- length(paid)
      if (length(paid)) toppred[y] <- preds[paid[which.max(bb[paid])]]
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
  if (!is.null(panel) && isTRUE(panel$ambiguous)) {
    # 0.3.0: balanced-panel symmetry -- the measures are the same under every
    # reading, so they ARE reportable; the entity/time axes are not.
    meas <- Reduce(intersect, lapply(panel$panels, `[[`, "measures"))
    axes <- lapply(panel$panels, function(p) c(cluster = p$cluster, time = p$time[1]))
    if (length(meas) >= 1L) {
      return(list(found = TRUE, role = "repeated_measures", column = meas,
                  detected_by = "mathematical", score = 70, max_score = 120, pct = 58.3,
                  panel_axes = axes,
                  reason = "balanced panel: the entity/time axis reading is value-symmetric and left to the user (see panel_axes)",
                  breakdown = list(list(name = "Long-format repeated measures (axis-ambiguous panel)",
                                        score = 70, max = 120,
                                        detail = sprintf("%d measure(s); %d admissible axis readings",
                                                         length(meas), length(axes))))))
    }
    panel <- NULL
  }
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
  # Measured-and-rejected alternatives: see the v0.3.0 development report (rolescry_v030/v030_dev_report.md, not distributed with the package).
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
  if (length(blocks) > 0) {
    # 0.3.0 sink veto: a block member MDL-explained by >= 2 MUTUALLY INDEPENDENT
    # other members is a dependency SINK (collider), not a repeated measurement --
    # genuine repeated measures share a latent factor, so their "causes" are never
    # mutually independent. Reuses the outcome machinery; no new constant. Probe D:
    # y = x1 + x2 + noise was swallowed as {x1, x2, y} and outcome_continuous never
    # fired. Members are value-fingerprint-ordered, so the sweep is equivariant.
    drop_sinks <- function(members) {
      keep <- members
      for (m in members) {
        others <- setdiff(keep, m)
        if (length(others) < 2L) next
        yv <- M[, m]; Xo <- M[, others, drop = FALSE]
        cc <- which(is.finite(yv) & rowSums(!is.finite(Xo)) == 0)
        if (length(cc) < 10L) next
        if (length(cc) > .MDL_NMAX) cc <- cc[.dcor_idx(length(cc), .MDL_NMAX, seed = 7L)]
        fpo <- vapply(others, function(c) .col_fingerprint(data[[c]]), character(1L))
        r <- .mdl_forward_gauss(yv[cc], Xo[cc, , drop = FALSE], fpo)
        if (r$k >= 2L && r$bits_saved > 0 &&
            .mdl_set_independent(Xo[cc, r$sel, drop = FALSE], fpo[r$sel]))
          keep <- setdiff(keep, m)
      }
      keep
    }
    blocks <- lapply(blocks, drop_sinks)
    blocks <- blocks[vapply(blocks, length, 0L) >= 3L]
    blocks <- blocks[!duplicated(vapply(blocks, paste, character(1L), collapse = ""))]
  }
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

# Scale-item column type: "binary" ({0,1} exactly), "ordinal" (integer-valued, >= 3 levels,
# discrete), else NA. Needs >= 20 observed values.
.scale_item_type <- function(x) {
  nz <- x[!is.na(x)]
  if (length(nz) < 20L || !is.numeric(x)) return(NA_character_)
  u <- sort(unique(nz), method = "radix")
  if (length(u) == 2L && isTRUE(all(abs(u - c(0, 1)) < 1e-9))) return("binary")
  if (all(nz == round(nz)) && length(u) >= 3L && .is_discrete(x)) return("ordinal")
  NA_character_
}

.detect_scale_items <- function(data, var_info) {
  # Ordinal items: the fixed candidate set must win the one-factor MDL code with a common
  # quantity (independent columns give negative bits). Binary batteries stay on the ported
  # scorer pending a Bernoulli one-latent code.
  none <- list(found = FALSE, role = "scale_items", column = NULL,
               detected_by = "mathematical", score = 0, max_score = NA_real_, pct = 0, breakdown = list())
  hit <- function(cols, score, nm, detail) list(
    found = TRUE, role = "scale_items", column = cols, detected_by = "mathematical",
    score = score, max_score = NA_real_, pct = max(0, round(100 * (1 - 2^(-score / 8)))),
    breakdown = list(list(name = nm, score = score, max = NA_real_, detail = detail)),
    n_items = length(cols))
  ord_cols <- var_info$column[var_info$user_type %in% c("Ordinal", "Continuous")]
  cand <- Filter(function(col) identical(.scale_item_type(data[[col]]), "ordinal"), ord_cols)
  if (length(cand) >= 3L) {
    f <- .mdl_factor_bits(as.matrix(data.frame(lapply(data[cand], as.numeric))))
    # same_quantity (communality > uniqueness) is a same-units / repeated-measures notion; scale items
    # are different questions sharing a latent, so it is reported but does not decide.
    if (f$bits > 0)
      return(hit(cand, f$bits, "one-factor MDL block",
                 sprintf("%d items; %.1f bits; same_quantity=%s", length(cand), f$bits, f$same_quantity)))
  }
  # D3 ANCESTOR PORT (complementary): a DICHOTOMOUS IRT battery -- >= 8 binary {0,1}
  # columns (value-only; the ordinal path above misses binaries). Only when the ordinal
  # path yielded no scale, so it never overrides the strong ordinal detector. STRUCTURAL
  # precision guard with NO tuned threshold: reject a one-hot PARTITION (every respondent
  # endorses exactly one item -> max row-sum <= 1), which is an encoding, not a scale; a
  # real battery has respondents endorsing several items (max row-sum > 1).
  {
    bin_cols <- character(0)
    for (col in var_info$column) {
      if (identical(.scale_item_type(data[[col]]), "binary")) bin_cols <- c(bin_cols, col)
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
  none
}

# 0.3.0: covariate is the RESIDUAL role -- a Continuous/Categorical/Binary column that no
# other role claimed. In 0.2.0 it listed every such column, including the survival pair,
# group and outcomes, so it carried no information beyond the column types.
# 0.3.0: calendar axis (time_index). Every Datetime-typed column is a calendar
# axis; the breakdown records whether it is monotone in stored order and regularly
# spaced. Needed for finance and longitudinal data, where the date column was
# previously typed "Other" and invisible to every detector.
.detect_time_index <- function(data, var_info) {
  tcols <- var_info$column[var_info$user_type == "Datetime"]
  if (length(tcols) == 0L) {
    return(list(found = FALSE, role = "time_index", column = NULL,
                detected_by = "mathematical", score = 0, max_score = 100, pct = 0,
                breakdown = list()))
  }
  bd <- lapply(tcols, function(col) {
    nv <- .datetime_num(data[[col]]); nv <- nv[is.finite(nv)]; dfs <- diff(nv)
    list(name = col, score = 80, max = 100,
         detail = sprintf("calendar values; monotone=%s regular_spacing=%s",
                          length(dfs) > 0 && all(dfs > 0),
                          length(dfs) > 1 && isTRUE(stats::sd(dfs) == 0)))
  })
  list(found = TRUE, role = "time_index", column = tcols,
       detected_by = "mathematical", score = 80, max_score = 100, pct = 80,
       breakdown = bd)
}

# Numeric view of a (possibly character-coded) calendar column.
.datetime_num <- function(v) {
  if (is.character(v) || is.factor(v)) {
    d <- tryCatch(suppressWarnings(as.Date(as.character(v), tryFormats = .DATE_FORMATS)),
                  error = function(e) as.Date(rep(NA, length(v))))
    as.numeric(d)
  } else suppressWarnings(as.numeric(v))
}

.detect_covariates <- function(data, var_info, claimed = character(0)) {
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  cat_cols <- var_info$column[var_info$user_type %in% c("Categorical", "Binary")]
  covs <- setdiff(c(num_cols, cat_cols), claimed)
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
  cols <- if (is.null(col)) character(0)
  else if (is.list(col)) unname(unlist(lapply(col, function(p) unlist(p))))
  else as.character(col)
  # 0.3.0: the axis candidates of an axis-ambiguous panel are claimed too -- a
  # column offered as an admissible entity/time reading is not a residual covariate.
  unique(c(cols, unname(unlist(r$panel_axes))))
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
    coupled_pair       = .detect_survival_components(data, var_info, "pair", name_bonus),
    time_index         = .detect_time_index(data, var_info),
    subject_id         = .detect_subject_id(data, var_info, name_bonus, claimed),
    outcome_continuous = .detect_outcome(data, var_info, "continuous", name_bonus, claimed),
    outcome_binary     = .detect_outcome(data, var_info, "binary", name_bonus, claimed),
    repeated_measures  = .detect_repeated_measures(data, var_info),
    scale_items        = .detect_scale_items(data, var_info),
    covariate          = .detect_covariates(data, var_info, claimed),
    list(found = FALSE, role = need, column = NULL, detected_by = "none",
         score = 0, max_score = 0, pct = 0, breakdown = list())
  )
}

# The canonical role-need set, in detection order. Structural roles (survival, subject id,
# repeated, scale) resolve BEFORE the outcomes so their columns are `claimed` (atomicity).
.ROLE_NEEDS <- c(
  "group_var", "paired_pairs", "agreement_pairs", "time_variable",
  "event_variable", "coupled_pair", "time_index", "subject_id",
  "repeated_measures", "scale_items", "outcome_continuous", "outcome_binary",
  "covariate"
)

# Roles whose detected columns BLOCK a column from also being an outcome (atomicity).
# Excludes covariate (the catch-all) and the outcome roles themselves.
.CLAIMING_ROLES <- c("group_var", "paired_pairs", "agreement_pairs", "time_variable",
                     "event_variable", "coupled_pair", "time_index", "subject_id",
                     "repeated_measures", "scale_items")

# Run every detector + per-column value classification + potential pairs.
run_all_detections <- function(data, var_info, name_bonus = NULL, verbose = FALSE) {
  roles <- list()
  for (need in .ROLE_NEEDS) {
    owners <- if (need == "covariate") names(roles) else intersect(names(roles), .CLAIMING_ROLES)
    # Paired/agreement members are not outcome candidates; when the paired scorer latches onto a
    # (cause, outcome) pair the continuous outcome stays silent (known limit).
    # Measured-and-rejected alternatives: see the v0.3.0 development report (rolescry_v030/v030_dev_report.md, not distributed with the package).
    claimed <- unique(unlist(lapply(owners, function(rn) .role_cols(roles[[rn]]))))
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
        mb <- .mdl_pair_bits(a, b)
        if (isTRUE(mb$bits > 0 && mb$agree)) {
          potential_pairs[[length(potential_pairs) + 1]] <- list(
            col1 = num_cols[i], col2 = num_cols[j],
            bits = mb$bits, shift_bits = .mdl_shift_bits(a, b)
          )
        }
      }
    }
  }
  list(roles = roles, value_types = value_types, potential_pairs = potential_pairs)
}
