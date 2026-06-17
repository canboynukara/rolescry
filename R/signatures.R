# Mathematical signature scorers. Ported verbatim (logic) from MDStatR
# engine/run_097_detect_roles.R, with stats:: qualification and ks.test wrapped
# in suppressWarnings() (the only baseline warning was the benign "p-value will
# be approximate in the presence of ties"). All scoring is value-based; no
# column names are read. No exports.

.score_paired_signature <- function(a, b) {
  bd <- list()
  n <- sum(stats::complete.cases(cbind(a, b)))
  if (n < 10) {
    return(list(score = 0, max = 110, pct = 0, breakdown = list()))
  }
  r <- tryCatch(stats::cor(a, b, use = "complete.obs"), error = function(e) NA)
  s1 <- if (!is.na(r) && r >= 0.3 && r <= 0.95) round(20 * min(1, (r - 0.1) / 0.6)) else 0L
  bd[[1]] <- list(name = "Correlation", score = s1, max = 20, detail = sprintf("r=%.2f", ifelse(is.na(r), 0, r)))

  cc <- stats::complete.cases(cbind(a, b))
  d <- mean(a[cc] - b[cc], na.rm = TRUE)
  pooled <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  ratio <- if (!is.na(pooled) && pooled > 0) abs(d) / pooled else Inf
  s2 <- if (is.finite(ratio) && ratio < 2) round(20 * (1 - ratio / 2)) else 0L
  bd[[2]] <- list(name = "Mean diff ratio", score = s2, max = 20, detail = sprintf("%.2f", ifelse(is.finite(ratio), ratio, Inf)))

  sd_a <- stats::sd(a, na.rm = TRUE)
  sd_b <- stats::sd(b, na.rm = TRUE)
  sdr <- if (!is.na(sd_a) && !is.na(sd_b) && sd_b > 1e-10) sd_a / sd_b else NA
  s3 <- if (!is.na(sdr) && sdr >= 0.3 && sdr <= 3) round(20 * (1 - abs(log(sdr)) / log(3))) else 0L
  bd[[3]] <- list(name = "SD similarity", score = s3, max = 20, detail = sprintf("ratio=%.2f", ifelse(is.na(sdr), 0, sdr)))

  rng_a <- range(a, na.rm = TRUE)
  rng_b <- range(b, na.rm = TRUE)
  ov <- max(0, min(rng_a[2], rng_b[2]) - max(rng_a[1], rng_b[1]))
  total_range <- max(rng_a[2], rng_b[2]) - min(rng_a[1], rng_b[1])
  ov_pct <- if (total_range > 0) ov / total_range else 0
  s4 <- if (ov_pct > 0.3) round(20 * min(1, ov_pct)) else 0L
  bd[[4]] <- list(name = "Range overlap", score = s4, max = 20, detail = sprintf("%.0f%%", ov_pct * 100))

  ks_p <- tryCatch(suppressWarnings(stats::ks.test(a[cc], b[cc])$p.value), error = function(e) 0)
  if (is.na(ks_p)) ks_p <- 0
  s5 <- if (ks_p > 0.05) round(15 * min(1, ks_p / 0.5)) else 0L
  bd[[5]] <- list(name = "Distribution shape", score = s5, max = 15, detail = sprintf("KS p=%.3f", ks_p))

  cpr <- n / max(length(a), 1)
  s6 <- if (cpr > 0.7) round(15 * min(1, (cpr - 0.5) / 0.5)) else 0L
  bd[[6]] <- list(name = "Complete pairs", score = s6, max = 15, detail = sprintf("%.0f%%", cpr * 100))

  total <- s1 + s2 + s3 + s4 + s5 + s6
  list(score = total, max = 110, pct = round(total / 110 * 100, 1), breakdown = bd)
}

.score_group_signature <- function(column) {
  bd <- list()
  if (is.null(column) || (!is.character(column) && !is.factor(column))) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  vals <- column
  vals <- vals[!is.na(vals)]
  if (length(vals) < 5) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  lvls <- unique(vals)
  n_levels <- length(lvls)

  s1 <- if (n_levels >= 2 && n_levels <= 10) {
    if (n_levels <= 5) 25L else round(25 * (1 - (n_levels - 5) / 10))
  } else {
    0L
  }
  bd[[1]] <- list(name = "Level count", score = s1, max = 25, detail = sprintf("%d levels", n_levels))

  counts <- table(vals)
  balance <- min(counts) / max(counts)
  s2 <- round(25 * min(1, balance / 0.5))
  bd[[2]] <- list(name = "Balance", score = s2, max = 25, detail = sprintf("min/max=%.2f", balance))

  props <- as.numeric(counts) / sum(counts)
  entropy <- -sum(props * log2(props + 1e-10))
  max_entropy <- log2(n_levels)
  norm_entropy <- if (max_entropy > 0) entropy / max_entropy else 0
  s3 <- if (norm_entropy > 0.5) round(20 * min(1, norm_entropy)) else 0L
  bd[[3]] <- list(name = "Entropy", score = s3, max = 20, detail = sprintf("%.2f", norm_entropy))

  is_nonnumeric <- !all(grepl("^-?[0-9.]+$", as.character(lvls)))
  s4 <- if (is_nonnumeric) {
    15L
  } else {
    if (n_levels <= 5) 10L else 0L
  }
  bd[[4]] <- list(name = "Non-numeric", score = s4, max = 15, detail = ifelse(is_nonnumeric, "TRUE", "FALSE"))

  na_pct <- mean(is.na(column)) * 100
  s5 <- if (na_pct < 5) 15L else if (na_pct < 20) round(15 * (1 - na_pct / 40)) else 0L
  bd[[5]] <- list(name = "Completeness", score = s5, max = 15, detail = sprintf("%.1f%% NA", na_pct))

  total <- s1 + s2 + s3 + s4 + s5
  list(
    score = total, max = 100, pct = round(total / 100 * 100, 1), breakdown = bd,
    details = list(n_levels = n_levels, level_counts = as.numeric(counts))
  )
}

.score_survival_signature <- function(time_col, event_col, data) {
  bd <- list()
  tv <- data[[time_col]]
  ev <- data[[event_col]]
  if (is.null(tv) || is.null(ev)) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  tv <- as.numeric(tv)
  ev <- as.numeric(as.character(ev))
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

.score_agreement_signature <- function(a, b) {
  bd <- list()
  n <- sum(stats::complete.cases(cbind(a, b)))
  if (n < 10) {
    return(list(score = 0, max = 110, pct = 0, breakdown = list()))
  }
  cc <- stats::complete.cases(cbind(a, b))
  ac <- a[cc]
  bc <- b[cc]

  r <- tryCatch(stats::cor(ac, bc), error = function(e) NA)
  s1 <- if (!is.na(r) && r >= 0.3 && r <= 0.85) round(15 * min(1, (r - 0.1) / 0.5)) else if (!is.na(r) && r > 0.85) round(15 * max(0, 1 - (r - 0.85) / 0.15)) else 0L
  bd[[1]] <- list(name = "Correlation", score = s1, max = 15, detail = sprintf("r=%.2f", ifelse(is.na(r), 0, r)))

  diffs <- ac - bc
  bias <- mean(diffs)
  rng <- max(ac, bc) - min(ac, bc)
  bias_ratio <- if (rng > 0) abs(bias) / rng else Inf
  s2 <- if (is.finite(bias_ratio) && bias_ratio < 0.3) round(25 * (1 - bias_ratio / 0.3)) else 0L
  bd[[2]] <- list(name = "BA bias/range", score = s2, max = 25, detail = sprintf("%.3f", ifelse(is.finite(bias_ratio), bias_ratio, Inf)))

  sd_a <- stats::sd(ac)
  sd_b <- stats::sd(bc)
  sdr <- if (!is.na(sd_b) && sd_b > 1e-10) sd_a / sd_b else NA
  s3 <- if (!is.na(sdr) && sdr >= 0.5 && sdr <= 2.0) round(20 * (1 - abs(log(sdr)) / log(2))) else 0L
  bd[[3]] <- list(name = "SD ratio", score = s3, max = 20, detail = sprintf("%.2f", ifelse(is.na(sdr), 0, sdr)))

  rng_a <- range(ac)
  rng_b <- range(bc)
  ov <- max(0, min(rng_a[2], rng_b[2]) - max(rng_a[1], rng_b[1]))
  tr <- max(rng_a[2], rng_b[2]) - min(rng_a[1], rng_b[1])
  ov_pct <- if (tr > 0) ov / tr else 0
  s4 <- if (ov_pct > 0.5) round(20 * min(1, ov_pct)) else 0L
  bd[[4]] <- list(name = "Range overlap", score = s4, max = 20, detail = sprintf("%.0f%%", ov_pct * 100))

  icc_val <- tryCatch(
    {
      mat <- cbind(ac, bc)
      k <- ncol(mat)
      ni <- nrow(mat)
      grand <- mean(mat)
      ss_r <- k * sum((rowMeans(mat) - grand)^2)
      ss_c <- ni * sum((colMeans(mat) - grand)^2)
      ss_t <- sum((mat - grand)^2)
      ss_e <- ss_t - ss_r - ss_c
      ms_r <- ss_r / (ni - 1)
      ms_e <- ss_e / ((ni - 1) * (k - 1))
      ms_c <- ss_c / (k - 1)
      (ms_r - ms_e) / (ms_r + (k - 1) * ms_e + k * (ms_c - ms_e) / ni)
    },
    error = function(e) NA
  )
  s5 <- if (!is.na(icc_val) && icc_val > 0.3) round(15 * min(1, icc_val)) else 0L
  bd[[5]] <- list(name = "ICC(2,1)", score = s5, max = 15, detail = sprintf("%.2f", ifelse(is.na(icc_val), 0, icc_val)))

  kurt_diff <- tryCatch(
    {
      k_a <- mean((ac - mean(ac))^4) / (stats::sd(ac)^4) - 3
      k_b <- mean((bc - mean(bc))^4) / (stats::sd(bc)^4) - 3
      abs(k_a - k_b)
    },
    error = function(e) Inf
  )
  s6 <- if (is.finite(kurt_diff) && kurt_diff < 2) round(15 * (1 - kurt_diff / 2)) else 0L
  bd[[6]] <- list(name = "Kurtosis similarity", score = s6, max = 15, detail = sprintf("diff=%.2f", ifelse(is.finite(kurt_diff), kurt_diff, Inf)))

  total <- s1 + s2 + s3 + s4 + s5 + s6
  list(score = total, max = 110, pct = round(total / 110 * 100, 1), breakdown = bd)
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
  cv_means <- if (mean(col_means) != 0) stats::sd(col_means) / abs(mean(col_means)) else Inf
  s3 <- if (is.finite(cv_means) && cv_means < 1) round(20 * (1 - cv_means)) else 0L
  bd[[3]] <- list(name = "Mean similarity", score = s3, max = 20, detail = sprintf("CV=%.2f", ifelse(is.finite(cv_means), cv_means, Inf)))

  col_sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
  cv_sds <- if (mean(col_sds) > 0) stats::sd(col_sds) / mean(col_sds) else Inf
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

.score_scale_items_signature <- function(col_list, data) {
  bd <- list()
  if (length(col_list) < 5) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }
  mat <- as.matrix(data[, col_list, drop = FALSE])
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  n <- nrow(mat)
  k <- ncol(mat)
  if (n < 20) {
    return(list(score = 0, max = 100, pct = 0, breakdown = list()))
  }

  s1 <- if (k >= 5) 20L else 0L
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
