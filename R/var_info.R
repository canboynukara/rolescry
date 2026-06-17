# Value-based column typing. Ported from MDStatR engine/build_rv_typed.R
# (.vb_build_var_info, F0887 "Data inspice, non nomen"): a column's type is
# decided by its VALUES. A small name hint for ID columns is value-vetoed (a
# binary or non-integer column under an id-ish name is typed by value, not ID).
#
# Returns data.frame(column, user_type) with user_type in
# {Binary, ID, Categorical, Continuous, Other}.
.build_var_info <- function(d) {
  .name_id <- function(nm) {
    grepl("(?i)(^|_)(id|participant|record(_id)?|case(_id)?)$", nm, perl = TRUE)
  }
  data.frame(
    column = names(d),
    user_type = vapply(seq_along(d), function(i) {
      x <- d[[i]]
      nm <- names(d)[i]
      name_id <- .name_id(nm)
      if (is.character(x) || is.factor(x)) {
        n <- length(x)
        xv <- x[!is.na(x)]
        uv <- length(unique(xv))
        if (uv == 2L) return("Binary")
        if (n >= 10L && uv == n) return("ID")
        if (name_id) return("ID")
        return("Categorical")
      }
      if (is.numeric(x)) {
        xv <- x[!is.na(x)]
        n <- length(x)
        uv <- length(unique(xv))
        if (uv == 2L) return("Binary")
        is_int <- length(xv) > 0L && all(abs(xv - round(xv)) < 1e-9)
        if (n >= 10L && uv == n && is_int) return("ID")
        if (name_id && is_int) return("ID")
        return("Continuous")
      }
      "Other"
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}
