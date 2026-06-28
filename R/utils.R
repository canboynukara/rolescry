# Internal utilities. No exports.

# NULL/length-0/NA-coalesce (canonical NA-aware variant, severed from MDStatR
# build_rv_typed.R / _guards.R where it was an injected global).
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (is.atomic(a) && length(a) == 1L && is.na(a)) return(b)
  a
}

# Locale-independent ASCII case fold. Base-only replacement for MDStatR's
# md_casefold_ascii (utils.R:1417), which relied on stringi or a process-wide
# Sys.setlocale() mutation (a CRAN --as-cran blocker). chartr over A-Z is
# locale-invariant, so this folds uppercase ASCII letters without touching the
# session locale or depending on stringi.
.casefold_ascii <- function(x, case = c("lower", "upper")) {
  case <- match.arg(case)
  if (length(x) == 0L) return(x)
  x <- as.character(x)
  if (case == "lower") {
    chartr("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz", x)
  } else {
    chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", x)
  }
}

# Capped (<= max_pts) name bonus for a single column name against a keyword set.
# Returns 0 when no keywords (name_bonus = NULL path) so detection stays purely
# mathematical by default. Keywords are matched case-insensitively as regex
# fragments on the ASCII-folded column name.
.name_bonus_pts <- function(colname, keywords, max_pts = 10) {
  if (is.null(keywords) || length(keywords) == 0L) return(0)
  cl <- .casefold_ascii(colname, "lower")
  hit <- any(vapply(keywords, function(k) grepl(k, cl), logical(1L)))
  if (isTRUE(hit)) max_pts else 0
}

# Conditional diagnostic emit (stderr via message), gated by verbose. Replaces
# every unconditional cat() in the MDStatR sources (a CRAN blocker).
.say <- function(verbose, ...) {
  if (isTRUE(verbose)) message(sprintf(...))
  invisible(NULL)
}

# Order-invariant value fingerprint of a column (v0.1.1, audit Finding C1).
# A deterministic string key computed from the column's VALUES (and fixed row
# positions, which never permute) and never from its column index. Used as the
# canonical tie-break among equal-scoring candidates so the winner is the same
# regardless of column order -> S_n (permutation) equivariance. Two distinct
# columns collide only if value-identical (a genuine data symmetry).
.col_fingerprint <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    lv <- sort(unique(x[!is.na(x)]), method = "radix")
    x <- as.integer(factor(x, levels = lv))            # value-sorted codes (name-free)
  }
  xn <- suppressWarnings(as.numeric(x))
  ok <- !is.na(xn)
  if (!any(ok)) return(sprintf("E|%d", length(xn)))
  v <- xn[ok]; w <- which(ok)                          # original row indices (stable under col perm)
  paste(length(xn), length(unique(v)),
        format(sum(v), digits = 17), format(sum(v * w), digits = 17),
        format(sum(v * v), digits = 17),
        format(min(v), digits = 17), format(max(v), digits = 17), sep = "|")
}

# Normalize a candidate survival-event column to {0,1} (v0.1.1, audit C3).
# v0.1.0 required the event coded EXACTLY {0,1}, rejecting the extremely common
# {1,2} (e.g. survival::lung), {0,1,2} (survival::pbc), and logical codings.
# Accepts logical, numeric 2-level ({0,1},{1,2},...), and 3-level integer status
# ({0,1,2}: 0=reference/censor, >0 => event). NA preserved.
# 2-level factor/character is deliberately NOT treated as a survival event:
# classification TARGETS are commonly factor-coded and would be mis-detected as
# survival events (false positives) -- clinical event indicators are numeric.
# Returns a {0,1} integer vector or NULL if not event-like.
.normalize_event <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (!is.numeric(x)) return(NULL)
  xn <- suppressWarnings(as.numeric(x))
  un <- sort(unique(xn[!is.na(xn)]), method = "radix")
  if (length(un) == 2L) return(as.integer(match(xn, un)) - 1L)        # min->0, max->1
  if (length(un) == 3L && all(abs(un - round(un)) < 1e-9))
    return(as.integer(xn > min(un)))                                  # {0,1,2}: 0=censor, >0=event
  NULL
}

# Canonical winner among candidate columns: maximize `scores`, break ties by
# smallest value-fingerprint. `cols` parallel to `scores`. Equivariant.
# Returns NA_character_ on an UNBREAKABLE tie -- when the top two candidates
# share both score AND value-fingerprint (value-identical columns). For such a
# genuine data symmetry no single pick can be both relabel-invariant and
# S_n-equivariant, so the honest, invariant answer is silence (caller treats
# NA as not-found). Distinct columns never collide barring identical values.
.canonical_pick <- function(cols, scores, data) {
  fps <- unname(vapply(cols, function(k) .col_fingerprint(data[[k]]), character(1L)))
  ord <- order(-scores, fps)
  if (length(ord) >= 2L) {
    i1 <- ord[1]; i2 <- ord[2]
    if (isTRUE(scores[i1] == scores[i2]) && isTRUE(fps[i1] == fps[i2])) return(NA_character_)
  }
  cols[ord[1]]
}
