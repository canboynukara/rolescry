# Mathematical signature scorers. Ported verbatim (logic) from MDStatR
# engine/run_097_detect_roles.R, with stats:: qualification and ks.test wrapped
# in suppressWarnings() (the only baseline warning was the benign "p-value will
# be approximate in the presence of ties"). All scoring is value-based; no
# column names are read. No exports.

.score_survival_signature <- function(time_col, event_col, data) {
  bd <- list()
  tv <- data[[time_col]]
  ev <- data[[event_col]]
  if (is.null(tv) || is.null(ev)) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  tv <- as.numeric(tv)
  ev <- .normalize_event(ev)                              # C3: normalize coding to {0,1}
  if (is.null(ev)) return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  n <- sum(stats::complete.cases(cbind(tv, ev)))
  if (n < 10) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  cc <- stats::complete.cases(cbind(tv, ev))
  tv_c <- tv[cc]
  ev_c <- ev[cc]

  uvals <- sort(unique(ev_c), method = "radix")
  s1 <- if (length(uvals) == 2 && all(uvals %in% c(0, 1))) 20L else 0L
  bd[[1]] <- list(name = "Binary event", score = s1, max = 20, detail = paste(uvals, collapse = ","))

  all_pos <- all(tv_c > 0, na.rm = TRUE)
  skew <- tryCatch(
    {
      m3 <- mean((tv_c - mean(tv_c))^3)
      m3 / (stats::sd(tv_c)^3)
    },
    error = function(e) 0
  )
  if (!is.finite(skew)) skew <- 0                       # v0.1.1: guard NaN/Inf (sd==0) -> no crash
  s2_pos <- if (all_pos) 10L else 0L
  s2_skew <- if (skew > 0.3) min(10L, round(10 * min(1, skew / 2))) else 0L
  s2 <- s2_pos + s2_skew
  bd[[2]] <- list(name = "Time positive+skew", score = s2, max = 20, detail = sprintf("pos=%s skew=%.2f", all_pos, skew))

  event_rate <- mean(ev_c == 1)
  s3 <- if (event_rate >= 0.05 && event_rate <= 0.95) {
    round(20 * (1 - abs(event_rate - 0.5) / 0.5))
  } else {
    0L
  }
  bd[[3]] <- list(name = "Event rate", score = s3, max = 20, detail = sprintf("%.1f%%", event_rate * 100))

  med_event <- stats::median(tv_c[ev_c == 1], na.rm = TRUE)
  med_censor <- stats::median(tv_c[ev_c == 0], na.rm = TRUE)
  s4 <- if (!is.na(med_event) && !is.na(med_censor) && med_event < med_censor) 20L else if (!is.na(med_event) && !is.na(med_censor)) 5L else 0L
  bd[[4]] <- list(
    name = "Event timing", score = s4, max = 20,
    detail = sprintf(
      "med_event=%.1f med_censor=%.1f", ifelse(is.na(med_event), 0, med_event),
      ifelse(is.na(med_censor), 0, med_censor)
    )
  )

  max_time <- max(tv_c, na.rm = TRUE)
  ties_at_max <- sum(abs(tv_c - max_time) < .Machine$double.eps^0.5 * max(1, abs(max_time)))
  s5 <- if (ties_at_max <= 3) 20L else if (ties_at_max <= n * 0.1) 10L else 0L
  bd[[5]] <- list(name = "Max-time ties", score = s5, max = 20, detail = sprintf("%d ties at max", ties_at_max))

  total <- s1 + s2 + s3 + s4 + s5
  list(
    score = total, max = 100, pct = round(total / 100 * 100, 1), breakdown = bd,
    event_col = event_col, time_col = time_col
  )
}

.score_repeated_signature <- function(col_list, data) {
  bd <- list()
  if (length(col_list) < 3) {
    return(list(score = 0, max = 120, pct = 0, breakdown = list()))
  }
  mat <- as.matrix(data[, col_list, drop = FALSE])
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  n <- nrow(mat)
  if (n < 10) {
    return(list(score = 0, max = 120, pct = 0, breakdown = list()))
  }
  k <- ncol(mat)

  s1 <- if (k >= 3) 20L else 0L
  bd[[1]] <- list(name = "Column count", score = s1, max = 20, detail = sprintf("%d columns", k))

  cormat <- tryCatch(stats::cor(mat, use = "complete.obs"), error = function(e) matrix(NA, k, k))
  r_vals <- cormat[upper.tri(cormat)]
  mean_r <- mean(r_vals, na.rm = TRUE)
  s2 <- if (!is.na(mean_r) && mean_r > 0.3) round(20 * min(1, (mean_r - 0.1) / 0.6)) else 0L
  bd[[2]] <- list(name = "Mean inter-r", score = s2, max = 20, detail = sprintf("r=%.2f", ifelse(is.na(mean_r), 0, mean_r)))

  col_means <- colMeans(mat, na.rm = TRUE)
  mcm <- mean(col_means)                                 # AUREX: guard NaN-in-if (mixed +/-Inf column)
  cv_means <- if (is.finite(mcm) && mcm != 0) stats::sd(col_means) / abs(mcm) else Inf
  s3 <- if (is.finite(cv_means) && cv_means < 1) round(20 * (1 - cv_means)) else 0L
  bd[[3]] <- list(name = "Mean similarity", score = s3, max = 20, detail = sprintf("CV=%.2f", ifelse(is.finite(cv_means), cv_means, Inf)))

  col_sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
  msd <- mean(col_sds)                                   # AUREX: guard NaN-in-if
  cv_sds <- if (is.finite(msd) && msd > 0) stats::sd(col_sds) / msd else Inf
  s4 <- if (is.finite(cv_sds) && cv_sds < 1) round(20 * (1 - cv_sds)) else 0L
  bd[[4]] <- list(name = "SD similarity", score = s4, max = 20, detail = sprintf("CV=%.2f", ifelse(is.finite(cv_sds), cv_sds, Inf)))

  trend_r <- tryCatch(stats::cor(col_means, seq_along(col_means), method = "spearman"), error = function(e) 0)
  s5 <- round(20 * abs(ifelse(is.na(trend_r), 0, trend_r)))
  bd[[5]] <- list(name = "Trend", score = s5, max = 20, detail = sprintf("rho=%.2f", ifelse(is.na(trend_r), 0, trend_r)))

  s6 <- 10L
  bd[[6]] <- list(name = "Subject ID", score = s6, max = 20, detail = "base score")

  total <- s1 + s2 + s3 + s4 + s5 + s6
  list(score = total, max = 120, pct = round(total / 120 * 100, 1), breakdown = bd)
}

# Used by the binary IRT path of .detect_scale_items only (pending a Bernoulli one-latent code).
.score_scale_items_signature <- function(col_list, data) {
  bd <- list()
  # 0.3.0: 3-4 item short batteries (e.g. multi-agency ratings) are scoreable; the
  # caller guards them with the pairwise 2-sigma dependence gate.
  if (length(col_list) < 3) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  mat <- as.matrix(data[, col_list, drop = FALSE])
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  n <- nrow(mat)
  k <- ncol(mat)
  if (n < 20) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }

  s1 <- if (k >= 5) 20L else if (k >= 3) 10L else 0L
  bd[[1]] <- list(name = "Item count", score = s1, max = 20, detail = sprintf("%d items", k))

  ranges <- apply(mat, 2, function(x) length(unique(x)))
  median_range <- stats::median(ranges)
  s2 <- if (median_range >= 3 && median_range <= 10) 20L else if (median_range >= 2) 10L else 0L
  bd[[2]] <- list(name = "Ordinal range", score = s2, max = 20, detail = sprintf("median %d levels", median_range))

  cormat <- tryCatch(stats::cor(mat, use = "complete.obs"), error = function(e) matrix(NA, k, k))
  r_vals <- cormat[upper.tri(cormat)]
  mean_r <- mean(r_vals, na.rm = TRUE)
  s3 <- if (!is.na(mean_r) && mean_r > 0.15 && mean_r < 0.8) round(20 * min(1, mean_r / 0.4)) else 0L
  bd[[3]] <- list(name = "Inter-item r", score = s3, max = 20, detail = sprintf("%.2f", ifelse(is.na(mean_r), 0, mean_r)))

  alpha_proxy <- tryCatch(
    {
      k_items <- ncol(mat)
      item_vars <- apply(mat, 2, stats::var, na.rm = TRUE)
      total_var <- stats::var(rowSums(mat, na.rm = TRUE), na.rm = TRUE)
      (k_items / (k_items - 1)) * (1 - sum(item_vars) / total_var)
    },
    error = function(e) NA
  )
  s4 <- if (!is.na(alpha_proxy) && alpha_proxy > 0.6) round(20 * min(1, alpha_proxy)) else 0L
  bd[[4]] <- list(name = "Alpha proxy", score = s4, max = 20, detail = sprintf("%.2f", ifelse(is.na(alpha_proxy), 0, alpha_proxy)))

  all_int <- all(mat == round(mat), na.rm = TRUE)
  s5 <- if (all_int) 20L else 5L
  bd[[5]] <- list(name = "Integer values", score = s5, max = 20, detail = ifelse(all_int, "TRUE", "FALSE"))

  total <- s1 + s2 + s3 + s4 + s5
  list(score = total, max = 100, pct = round(total / 100 * 100, 1), breakdown = bd, n_items = k)
}
