.tidy_role <- function(r) {
  col <- r$column
  cols <- if (is.null(col)) {
    character(0)
  } else if (is.list(col)) {
    unname(unlist(lapply(col, function(p) unlist(p))))
  } else {
    as.character(col)
  }
  list(
    found = isTRUE(r$found),
    columns = cols,
    score = r$score %||% NA_real_,
    max_score = r$max_score %||% NA_real_,
    pct = r$pct %||% NA_real_,
    detected_by = r$detected_by %||% NA_character_,
    components = r$breakdown %||% list()
  )
}

#' Detect variable roles by data signature, not by name
#'
#' Inspects an already-loaded data frame and assigns each column (or group of
#' columns) to statistical roles -- group variable, continuous/binary outcome,
#' survival time and event, paired and agreement measurement pairs, repeated
#' measures, scale items, subject id, and covariates -- using only the data's
#' information-theoretic signature (Shannon entropy, distributional shape,
#' inter-column structure) and never the column names. Renaming columns to
#' \code{col_1, col_2, ...} does not change the result (the name-blindness, or
#' "turnusol", invariant).
#'
#' Detection is purely mathematical by default (\code{name_bonus = NULL}). When
#' a keyword dictionary is supplied via \code{name_bonus}, column names act only
#' as a small, capped tie-breaker (at most a +10 point nudge, i.e. <= 10 percent,
#' applied to candidate selection for the group, outcome, subject-id and
#' survival roles) -- the reported confidence stays the mathematical signature.
#' See [rolescry_default_name_bonus()] for a ready-made dictionary.
#'
#' @param data A \code{data.frame} (already loaded; for header-aware loading see
#'   [read_data()]).
#' @param name_bonus Optional named list mapping role keys to character vectors
#'   of case-insensitive keyword regex fragments, e.g.
#'   \code{list(group_var = c("treat", "arm"), outcome_binary = c("death"))}.
#'   Recognised keys: \code{group_var}, \code{outcome_continuous},
#'   \code{outcome_binary}, \code{subject_id}, \code{time_variable},
#'   \code{event_variable}. \code{NULL} (default) = pure signature detection.
#' @param verbose Logical; if \code{TRUE}, emit per-role progress via
#'   \code{message()}. Default \code{FALSE} (silent).
#'
#' @return An S3 object of class \code{"role_detection"}: a list with
#'   \describe{
#'     \item{var_info}{data.frame(column, type) -- value-based column typing.}
#'     \item{roles}{named list; each entry has \code{found}, \code{columns},
#'       \code{score}, \code{max_score}, \code{pct}, \code{detected_by} and a
#'       \code{components} score breakdown.}
#'     \item{value_types}{named character vector of per-column value-type labels.}
#'     \item{potential_pairs}{list of candidate continuous column pairs with
#'       paired and agreement scores.}
#'     \item{n_obs, n_var}{dataset dimensions.}
#'   }
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   arm  = rep(c(0, 1), each = 50),
#'   pre  = rnorm(100, 10, 2),
#'   post = rnorm(100, 11, 2),
#'   resp = rbinom(100, 1, 0.4)
#' )
#' res <- detect_roles(d)
#' res
#' res$roles$group_var$columns
#'
#' @seealso [read_data()], [compute_nmi()], [rolescry_default_name_bonus()]
#' @export
detect_roles <- function(data, name_bonus = NULL, verbose = FALSE) {
  if (!is.data.frame(data)) stop("detect_roles(): `data` must be a data.frame.")
  if (ncol(data) < 1L) stop("detect_roles(): `data` has no columns.")
  if (!is.null(name_bonus) && !is.list(name_bonus)) {
    stop("detect_roles(): `name_bonus` must be NULL or a named list of keyword vectors.")
  }
  vi <- .build_var_info(data)
  res <- run_all_detections(data, vi, name_bonus = name_bonus, verbose = verbose)
  out <- list(
    var_info = stats::setNames(
      data.frame(column = vi$column, type = vi$user_type, stringsAsFactors = FALSE),
      c("column", "type")
    ),
    roles = lapply(res$roles, .tidy_role),
    value_types = res$value_types,
    potential_pairs = res$potential_pairs,
    n_obs = nrow(data),
    n_var = ncol(data)
  )
  class(out) <- "role_detection"
  out
}

#' @param x A \code{role_detection} object.
#' @param ... Ignored.
#' @rdname detect_roles
#' @export
print.role_detection <- function(x, ...) {
  cat(sprintf("<role_detection> %d observations x %d variables\n", x$n_obs, x$n_var))
  found <- Filter(function(r) isTRUE(r$found), x$roles)
  if (length(found) == 0) {
    cat("  no roles detected\n")
  } else {
    for (nm in names(found)) {
      r <- found[[nm]]
      cat(sprintf(
        "  %-18s %-28s pct=%.1f%s\n", nm,
        paste(r$columns, collapse = ", "),
        r$pct,
        if (!is.na(r$detected_by) && grepl("name", r$detected_by)) " (+name)" else ""
      ))
    }
  }
  invisible(x)
}

#' @param object A \code{role_detection} object.
#' @rdname detect_roles
#' @export
summary.role_detection <- function(object, ...) {
  data.frame(
    role = names(object$roles),
    found = vapply(object$roles, function(r) isTRUE(r$found), logical(1L)),
    columns = vapply(object$roles, function(r) paste(r$columns, collapse = ","), character(1L)),
    pct = vapply(object$roles, function(r) as.numeric(r$pct), numeric(1L)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
