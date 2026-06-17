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
