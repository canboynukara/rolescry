# Universal mathematical value-type classifier. Ported verbatim (logic) from
# MDStatR engine/categorical_detector.R, with: stats:: qualification, the
# global-environment source-once sentinel removed, and ASCII-only comments.
# No dictionaries, no language detection. Works identically on col_N and named
# columns because it only ever sees the VECTOR, never the column name.
#
# Return shape: list(type, ...). Possible types: EMPTY, CONSTANT,
# BINARY_CATEGORICAL, BINARY_NUMERIC_ENCODED, ORDINAL_NUMERIC,
# CONTINUOUS_NUMERIC, CATEGORICAL_SURVEY_RESPONSE, NOMINAL_WITH_TYPOS,
# FREE_TEXT, NOMINAL_CATEGORICAL, CATEGORICAL_LIKELY, AMBIGUOUS_TEXT.
classify_value_type <- function(x) {
  x_nm <- x[!is.na(x) & x != ""]
  n <- length(x_nm)
  if (n == 0) {
    return(list(type = "EMPTY"))
  }
  u <- unique(x_nm)
  m <- length(u)
  if (m == 1) {
    return(list(type = "CONSTANT", n_unique = 1L))
  }
  if (m == 2 && !is.numeric(x)) {
    return(list(
      type = "BINARY_CATEGORICAL", n_unique = 2L,
      levels = sort(u, method = "radix")
    ))
  }
  if (is.numeric(x)) {
    finite_int_range <- is.finite(x_nm) & abs(x_nm) <= .Machine$integer.max
    int_ok <- isTRUE(all(finite_int_range) &&
      all(x_nm == suppressWarnings(as.integer(x_nm))))
    if (int_ok && m == 2) {
      return(list(
        type = "BINARY_NUMERIC_ENCODED", n_unique = 2L,
        levels = sort(u, method = "radix")
      ))
    }
    if (int_ok && m <= 9 && m >= 3) {
      return(list(type = "ORDINAL_NUMERIC", n_unique = m, levels = sort(u, method = "radix")))
    }
    return(list(type = "CONTINUOUS_NUMERIC", n_unique = m))
  }
  freq <- as.numeric(table(x_nm)) / n
  freq_pos <- freq[freq > 0]
  H <- -sum(freq_pos * log2(freq_pos))
  H_max <- log2(m)
  H_norm <- H / H_max
  uniq_ratio <- m / n
  med_len <- stats::median(nchar(as.character(x_nm)))
  freq_sorted <- sort(table(x_nm), decreasing = TRUE, method = "radix")
  k_for_80 <- which(cumsum(freq_sorted) / n >= 0.80)[1]
  if (k_for_80 <= 20 && uniq_ratio <= 0.05) {
    return(list(
      type = "CATEGORICAL_SURVEY_RESPONSE",
      top_k = as.integer(k_for_80),
      top_levels = names(freq_sorted)[1:k_for_80],
      tail_fraction = 1 - sum(freq_sorted[1:k_for_80]) / n,
      H_norm = H_norm
    ))
  }
  if (m > 7 && m <= 50 && requireNamespace("stringdist", quietly = TRUE)) {
    d <- stringdist::stringdistmatrix(u, u, method = "lv")
    cluster <- stats::cutree(stats::hclust(stats::as.dist(d), method = "single"), h = 2)
    n_cluster <- length(unique(cluster))
    if (n_cluster <= 7 && uniq_ratio <= 0.20) {
      return(list(
        type = "NOMINAL_WITH_TYPOS",
        n_raw = m, n_canonical = n_cluster, H_norm = H_norm
      ))
    }
  }
  if (uniq_ratio > 0.95) {
    return(list(type = "FREE_TEXT", reason = "uniq_ratio>0.95"))
  }
  if (med_len > 50) {
    return(list(type = "FREE_TEXT", reason = "med_len>50"))
  }
  if (m <= 20 && uniq_ratio <= 0.05) {
    return(list(
      type = "NOMINAL_CATEGORICAL", n_unique = m,
      levels = sort(u, method = "radix"), H_norm = H_norm
    ))
  }
  if (m <= 7 && uniq_ratio <= 0.20) {
    return(list(
      type = "CATEGORICAL_LIKELY", n_unique = m,
      levels = sort(u, method = "radix"), H_norm = H_norm
    ))
  }
  list(
    type = "AMBIGUOUS_TEXT", n_unique = m, uniq_ratio = uniq_ratio,
    H_norm = H_norm, med_len = med_len
  )
}
# (v0.1.1 AUREX) Removed dead `score_gap_ok` helper: defined but never called by
# any detector (the <=10% name-bonus cap is enforced directly in .name_bonus_pts
# via max_pts=10). Dead-code removal; no behavioral change.
