# Header-row detection. Ported from MDStatR engine/_pre_universal_loader.R
# (detect_header_row + extract_header + the .shannon_entropy / .char_class /
# .uniqueness / .length_score scorers), with: the unconditional cat() removed
# (verbose-gated), ASCII-only source, base+stats only (no file-format deps),
# and the weight config inlined as a constant. This is the pure scorer; file
# reading lives in read_data().

.HEADER_CFG <- list(
  W_ALPHA = 25, W_NONNUMERIC = 20, W_UNIQUE = 15, W_ENTROPY = 15,
  W_LENGTH = 10, W_TRANSITION = 15, W_COMPLETENESS = 15,
  MAX_SCAN_ROWS = 30, MIN_FILL_RATIO = 0.15
)

.shannon_entropy <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nchar(x) > 0]
  if (length(x) == 0) return(0)
  freq <- table(x)
  p <- as.numeric(freq) / sum(freq)
  p <- p[p > 0]
  if (length(p) <= 1) return(0)
  H <- -sum(p * log2(p))
  H_max <- log2(length(p))
  if (H_max == 0) return(0)
  H / H_max
}

.is_alpha_cell <- function(x) {
  # A cell is "alphabetic" when as.numeric() fails AND it contains any Unicode
  # letter. \p{L} (perl) covers Latin/accented/Cyrillic/etc. with ASCII source.
  x <- as.character(x)
  !is.na(x) & nchar(x) > 0 &
    suppressWarnings(is.na(as.numeric(x))) &
    grepl("\\p{L}", x, perl = TRUE)
}

.char_class_distribution <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nchar(x) > 0]
  if (length(x) == 0) return(list(alpha = 0, numeric = 0))
  n <- length(x)
  alpha_count <- sum(.is_alpha_cell(x))
  numeric_count <- sum(suppressWarnings(!is.na(as.numeric(x))))
  list(alpha = alpha_count / n, numeric = numeric_count / n)
}

.uniqueness <- function(x) {
  x <- x[!is.na(x) & nchar(as.character(x)) > 0]
  if (length(x) == 0) return(0)
  length(unique(x)) / length(x)
}

.length_score <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nchar(x) > 0]
  if (length(x) == 0) return(0)
  med_len <- stats::median(nchar(x))
  raw <- 1 / (1 + exp(-0.25 * (med_len - 5)))
  if (med_len > 50) raw <- raw * 0.3
  raw
}

.detect_header_row <- function(raw, verbose = FALSE) {
  if (is.null(raw) || nrow(raw) == 0 || ncol(raw) == 0) {
    return(list(header_row = 1, score = 0, all_scores = numeric(0), method = "fallback_empty"))
  }
  n_rows <- min(nrow(raw), .HEADER_CFG$MAX_SCAN_ROWS)
  n_cols <- ncol(raw)
  scores <- rep(-Inf, n_rows)
  for (i in seq_len(n_rows)) {
    row_vals <- as.character(unlist(raw[i, ]))
    non_na_mask <- !is.na(row_vals) & nchar(row_vals) > 0
    n_non_na <- sum(non_na_mask)
    if (n_non_na < max(3, n_cols * .HEADER_CFG$MIN_FILL_RATIO)) next
    vals <- row_vals[non_na_mask]
    cc <- .char_class_distribution(vals)
    s1 <- cc$alpha * .HEADER_CFG$W_ALPHA
    s2 <- (1 - cc$numeric) * .HEADER_CFG$W_NONNUMERIC
    s3 <- .uniqueness(vals) * .HEADER_CFG$W_UNIQUE
    s4 <- .shannon_entropy(vals) * .HEADER_CFG$W_ENTROPY
    s5 <- .length_score(vals) * .HEADER_CFG$W_LENGTH
    s6 <- 0
    if (i < nrow(raw)) {
      next_vals <- as.character(unlist(raw[i + 1, ]))
      next_vals <- next_vals[!is.na(next_vals) & nchar(next_vals) > 0]
      if (length(next_vals) > 0) {
        cc_next <- .char_class_distribution(next_vals)
        delta <- cc$alpha - cc_next$alpha
        if (!is.na(delta) && delta > 0) s6 <- min(delta, 1) * .HEADER_CFG$W_TRANSITION
      }
    }
    s7 <- (n_non_na / n_cols) * .HEADER_CFG$W_COMPLETENESS
    scores[i] <- s1 + s2 + s3 + s4 + s5 + s6 + s7
  }
  valid <- which(is.finite(scores))
  best_row <- if (length(valid) == 0) 1 else valid[which.max(scores[valid])]
  .say(verbose, "[rolescry] header row %d (score=%.1f)", best_row,
       if (is.finite(scores[best_row])) scores[best_row] else 0)
  list(
    header_row = best_row,
    score = if (is.finite(scores[best_row])) scores[best_row] else 0,
    all_scores = scores, method = "7signal"
  )
}

.extract_header <- function(raw, best_row) {
  n_cols <- ncol(raw)
  header_vals <- as.character(unlist(raw[best_row, ]))
  for (j in seq_len(n_cols)) {
    if (is.na(header_vals[j]) || nchar(trimws(header_vals[j])) == 0) {
      if (best_row > 1) {
        for (r in seq.int(best_row - 1L, 1L, by = -1L)) {
          parent_val <- as.character(raw[r, j])
          if (!is.na(parent_val) && nchar(trimws(parent_val)) > 0) {
            header_vals[j] <- trimws(parent_val)
            break
          }
        }
      }
    }
  }
  for (j in seq_len(n_cols)) {
    if (is.na(header_vals[j]) || nchar(trimws(header_vals[j])) == 0) {
      header_vals[j] <- paste0("col_", j)
    }
  }
  make.unique(trimws(header_vals), sep = "_")
}

#' Detect the header row of a raw, unparsed table
#'
#' Given a raw table read with \emph{no} header (every cell character), scores
#' each of the first rows with a 7-signal weighted heuristic (alphabetic ratio,
#' non-numeric ratio, uniqueness, normalized Shannon entropy, median string
#' length, alpha-vs-next-row transition, fill completeness) and returns the
#' most header-like row plus repaired, unique column names. Empty cells
#' are upward-filled (merged-cell repair) and any still-empty name becomes
#' \code{col_<j>}. Base + stats only; no file-format dependencies.
#'
#' @param raw A data.frame or matrix of the raw sheet, read with
#'   \code{header = FALSE} so the header row appears as data.
#' @param verbose Logical; emit the chosen row via \code{message()} if TRUE.
#'
#' @return A list with \code{header_row} (integer), \code{score} (numeric),
#'   \code{names} (repaired character vector, \code{length == ncol(raw)}),
#'   and \code{all_scores}.
#'
#' @examples
#' raw <- data.frame(
#'   V1 = c("age", "34", "51"),
#'   V2 = c("sex", "M", "F"),
#'   V3 = c("score", "8.1", "7.4"),
#'   stringsAsFactors = FALSE
#' )
#' detect_header(raw)$names
#'
#' @seealso [read_data()]
#' @export
detect_header <- function(raw, verbose = FALSE) {
  hr <- .detect_header_row(raw, verbose = verbose)
  nm <- .extract_header(raw, hr$header_row)
  list(header_row = hr$header_row, score = hr$score, names = nm, all_scores = hr$all_scores)
}
