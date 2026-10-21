# Value-based column typing. Ported from MDStatR engine/build_rv_typed.R
# (.vb_build_var_info, F0887 "Data inspice, non nomen"): a column's type is
# decided by its VALUES.
#
# v0.1.1 (audit Finding A): the name-based ID hint was REMOVED. In v0.1.0 a
# column whose NAME matched an id-ish regex was typed "ID" even under
# name_bonus=NULL (var_info.R:24,34), leaking the column name and breaking the
# TURNUSOL invariant for signature-weak (non-unique) ids. ID typing is now
# VALUE-ONLY (uniqueness + integer/character signature); any name hint flows
# solely through the name_bonus channel in the detectors, which is a no-op when
# name_bonus=NULL. Accepted consequence: a non-unique column under an id-ish
# name is no longer typed ID (it is Continuous/Categorical by value) -- correct,
# not a regression.
#
# Returns data.frame(column, user_type) with user_type in
# {Binary, ID, Categorical, Continuous, Other}.
# Accepted character date formats (value-only; TR day-first included).
.DATE_FORMATS <- c("%Y-%m-%d", "%Y/%m/%d", "%d.%m.%Y")

# All-values-parse date check for character columns: conservative (every non-NA
# value must parse), so a stray numeric or free-text column never types Datetime.
.parses_as_date <- function(xv) {
  if (length(xv) < 10L) return(FALSE)
  dt <- tryCatch(suppressWarnings(as.Date(xv, tryFormats = .DATE_FORMATS)),
                 error = function(e) NULL)
  !is.null(dt) && !anyNA(dt)
}

.build_var_info <- function(d) {
  data.frame(
    column = names(d),
    user_type = vapply(seq_along(d), function(i) {
      x <- d[[i]]
      # 0.3.0: calendar columns. The CLASS and the value FORMAT are part of the
      # value, not the name, so this preserves name-blindness.
      if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) return("Datetime")
      # 0.3.0: a logical column is a two-valued event/indicator coding (TRUE/FALSE),
      # not "Other". Typing it Other excluded every logical event column from the
      # survival and outcome candidate sets (pure bug, found by probe E).
      if (is.logical(x)) {
        uv <- length(unique(x[!is.na(x)]))
        return(if (uv == 2L) "Binary" else "Other")
      }
      if (is.character(x) || is.factor(x)) {
        n <- length(x)
        xv <- as.character(x[!is.na(x)])
        if (.parses_as_date(xv)) return("Datetime")
        uv <- length(unique(xv))
        if (uv == 2L) return("Binary")
        if (n >= 10L && uv == n) return("ID")
        return("Categorical")
      }
      if (is.numeric(x)) {
        xv <- x[!is.na(x)]
        n <- length(x)
        uv <- length(unique(xv))
        if (uv == 2L) return("Binary")
        is_int <- length(xv) > 0L && all(abs(xv - round(xv)) < 1e-9)
        if (n >= 10L && uv == n && is_int) return("ID")
        return("Continuous")
      }
      "Other"
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}
