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
  candidates <- var_info$column[var_info$user_type %in% c("Categorical", "Binary")]
  best_col <- NULL
  best_gs <- NULL
  best_rank <- -Inf
  for (col in candidates) {
    gs <- .score_group_signature(data[[col]])
    if (gs$pct >= 30) {
      rank <- gs$pct + .name_bonus_pts(col, name_bonus$group_var)
      if (rank > best_rank) {
        best_rank <- rank
        best_col <- col
        best_gs <- gs
      }
    }
  }
  if (!is.null(best_col)) {
    gs <- best_gs
    return(list(
      found = TRUE, role = "group_var", column = best_col,
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
      sc <- scorer(data[[num_cols[i]]], data[[num_cols[j]]])
      if (sc$pct >= 30) {
        pairs[[length(pairs) + 1]] <- list(col1 = num_cols[i], col2 = num_cols[j], score = sc)
      }
    }
  }
  if (length(pairs) > 0) {
    pairs <- pairs[order(-vapply(pairs, function(p) p$score$pct, numeric(1L)))]
    top <- pairs[[1]]
    return(list(
      found = TRUE, role = role,
      column = list(list(col1 = top$col1, col2 = top$col2)),
      all_pairs = pairs, detected_by = "mathematical",
      score = top$score$score, max_score = top$score$max, pct = top$score$pct,
      breakdown = top$score$breakdown
    ))
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
  best_time <- NULL
  best_event <- NULL
  best_ss <- NULL
  best_rank <- -Inf
  for (ev in binary_cols) {
    ev_vals <- suppressWarnings(as.numeric(as.character(data[[ev]])))
    if (any(is.na(ev_vals))) next
    if (length(unique(ev_vals)) != 2) next
    if (!all(sort(unique(ev_vals), method = "radix") == c(0, 1))) next
    for (tm in numeric_cols) {
      tm_vals <- suppressWarnings(as.numeric(data[[tm]]))
      if (any(tm_vals <= 0, na.rm = TRUE)) next
      ss <- .score_survival_signature(tm, ev, data)
      if (ss$pct >= 30) {
        rank <- ss$score + .name_bonus_pts(tm, name_bonus$time_variable) +
          .name_bonus_pts(ev, name_bonus$event_variable)
        if (rank > best_rank) {
          best_rank <- rank
          best_time <- tm
          best_event <- ev
          best_ss <- ss
        }
      }
    }
  }
  if (!is.null(best_ss)) {
    ss <- best_ss
    col <- if (component == "time") best_time else best_event
    return(list(
      found = TRUE, role = paste0(component, "_variable"), column = col,
      detected_by = "mathematical", score = ss$score, max_score = ss$max,
      pct = ss$pct, breakdown = ss$breakdown,
      event_col = best_event, time_col = best_time
    ))
  }
  list(
    found = FALSE, role = paste0(component, "_variable"), column = NULL,
    detected_by = "mathematical", score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

.detect_subject_id <- function(data, var_info, name_bonus = NULL) {
  id_cols <- var_info$column[var_info$user_type == "ID"]
  if (length(id_cols) > 0) {
    return(list(
      found = TRUE, role = "subject_id", column = id_cols[1],
      detected_by = "mathematical", score = 90, max_score = 100, pct = 90,
      breakdown = list(list(name = "ID type", score = 90, max = 100, detail = "tagged as ID")),
      needs_confirmation = FALSE
    ))
  }
  id_like_cols <- character(0)
  for (col in var_info$column) {
    vals <- data[[col]]
    if (!(length(unique(vals)) == nrow(data) && !any(is.na(vals)))) next
    id_like <- is.character(vals) || is.factor(vals) ||
      (is.numeric(vals) && isTRUE(all(vals == round(vals))))
    if (id_like) id_like_cols <- c(id_like_cols, col)
  }
  if (length(id_like_cols) > 0) {
    pick <- id_like_cols[1]
    for (col in id_like_cols) {
      if (.name_bonus_pts(col, name_bonus$subject_id) > 0) { pick <- col; break }
    }
    return(list(
      found = TRUE, role = "subject_id", column = pick,
      detected_by = "mathematical", score = 70, max_score = 100, pct = 70,
      breakdown = list(list(
        name = "Unique ID-like values", score = 70, max = 100,
        detail = "all unique, integer/character (not numeric-continuous)"
      )),
      needs_confirmation = TRUE
    ))
  }
  list(
    found = FALSE, role = "subject_id", column = NULL, detected_by = "mathematical",
    score = 0, max_score = 100, pct = 0, breakdown = list()
  )
}

.detect_outcome <- function(data, var_info, type, name_bonus = NULL) {
  if (type == "continuous") {
    cols <- var_info$column[var_info$user_type == "Continuous"]
  } else {
    cols <- var_info$column[var_info$user_type %in% c("Binary", "Categorical")]
    cols <- cols[vapply(cols, function(col_) length(unique(data[[col_]][!is.na(data[[col_]])])) == 2, logical(1L))]
  }
  if (length(cols) > 0) {
    pick <- cols[1]
    kw <- name_bonus[[paste0("outcome_", type)]]
    for (col in cols) {
      if (.name_bonus_pts(col, kw) > 0) { pick <- col; break }
    }
    return(list(
      found = TRUE, role = paste0("outcome_", type), column = pick,
      detected_by = "mathematical", score = 60, max_score = 100, pct = 60,
      breakdown = list(), needs_confirmation = TRUE
    ))
  }
  list(found = FALSE, role = paste0("outcome_", type), column = NULL,
       detected_by = "mathematical", score = 0, max_score = 100, pct = 0, breakdown = list())
}

.detect_repeated_measures <- function(data, var_info) {
  num_cols <- var_info$column[var_info$user_type == "Continuous"]
  if (length(num_cols) < 3) {
    return(list(found = FALSE, role = "repeated_measures", column = NULL,
                detected_by = "mathematical", score = 0, max_score = 120, pct = 0, breakdown = list()))
  }
  best_cols <- NULL
  best_rs <- NULL
  best_score <- 0
  for (start in seq_len(length(num_cols) - 2L)) {
    for (end_idx in (start + 2L):min(length(num_cols), start + 10L)) {
      cols <- num_cols[start:end_idx]
      rs <- .score_repeated_signature(cols, data)
      if (rs$pct >= 30 && rs$score > best_score) {
        best_cols <- cols
        best_rs <- rs
        best_score <- rs$score
      }
    }
  }
  if (best_score > 0) {
    rs <- best_rs
    return(list(found = TRUE, role = "repeated_measures", column = best_cols,
                detected_by = "mathematical", score = rs$score, max_score = rs$max,
                pct = rs$pct, breakdown = rs$breakdown))
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

# Dispatch one role need to its detector.
detect_variable_role <- function(data, var_info, need, name_bonus = NULL) {
  switch(need,
    group_var          = .detect_group_var(data, var_info, name_bonus),
    paired_pairs       = .detect_paired_pairs(data, var_info),
    agreement_pairs    = .detect_agreement_pairs(data, var_info),
    time_variable      = .detect_survival_components(data, var_info, "time", name_bonus),
    event_variable     = .detect_survival_components(data, var_info, "event", name_bonus),
    subject_id         = .detect_subject_id(data, var_info, name_bonus),
    outcome_continuous = .detect_outcome(data, var_info, "continuous", name_bonus),
    outcome_binary     = .detect_outcome(data, var_info, "binary", name_bonus),
    repeated_measures  = .detect_repeated_measures(data, var_info),
    scale_items        = .detect_scale_items(data, var_info),
    covariate          = .detect_covariates(data, var_info),
    list(found = FALSE, role = need, column = NULL, detected_by = "none",
         score = 0, max_score = 0, pct = 0, breakdown = list())
  )
}

# The canonical role-need set, in detection order.
.ROLE_NEEDS <- c(
  "group_var", "paired_pairs", "agreement_pairs", "time_variable",
  "event_variable", "subject_id", "outcome_continuous", "outcome_binary",
  "repeated_measures", "scale_items", "covariate"
)

# Run every detector + per-column value classification + potential pairs.
run_all_detections <- function(data, var_info, name_bonus = NULL, verbose = FALSE) {
  roles <- list()
  for (need in .ROLE_NEEDS) {
    roles[[need]] <- tryCatch(
      detect_variable_role(data, var_info, need, name_bonus),
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
