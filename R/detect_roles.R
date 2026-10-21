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
    components = r$breakdown %||% list(),
    confidence = r$confidence %||% NA_real_,    # v0.1.1: relational-outcome calibrated confidence (B1)
    reason = r$reason %||% NA_character_,        # v0.1.1: identifiability-silence explanation
    # 0.3.0 coupled_pair extras (NULL for every other role): which member is
    # which, and the two admissible readings the user chooses between.
    continuous_col = r$continuous_col %||% NULL,
    binary_col = r$binary_col %||% NULL,
    interpretations = r$interpretations %||% NULL,
    # 0.3.0 axis-ambiguous balanced panel: admissible (entity, time) readings.
    panel_axes = r$panel_axes %||% NULL
  )
}

#' Detect variable roles by data signature, not by name
#'
#' Inspects an already-loaded data frame and assigns each column (or group of
#' columns) to statistical roles -- group variable, continuous/binary outcome,
#' survival time and event, paired and agreement measurement pairs, repeated
#' measures, scale items, subject id, and covariates -- using only the data
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
#' A column can hold more than one role (e.g. a measurement in both
#' \code{paired_pairs} and \code{agreement_pairs}). The \code{covariate} role is
#' residual: it lists the continuous, categorical and binary columns that no other
#' role claimed. A role whose evidence does not single out a column is reported as
#' not found rather than guessed.
#'
#' The \code{coupled_pair} role (0.3.0) reports a dependent pair of one binary and
#' one positive continuous column whose reading is ambiguous: it may be a survival
#' (time, event) pair or a (cause, binary outcome) pair, and the data alone cannot
#' say which. The survival roles are assigned instead only when censoring geometry
#' singles the survival reading out. The entry carries \code{continuous_col},
#' \code{binary_col} and the two \code{interpretations}; the choice is the user's.
#'
#' The \code{time_index} role (0.3.0) reports calendar columns: Date/POSIXct
#' classes, or character columns whose every value parses as a date. The class
#' and format are properties of the values, so name-blindness is preserved. On a
#' balanced long panel whose entity/time axes are value-symmetric, the axes are
#' reported as ambiguous (\code{repeated_measures$panel_axes}) unless serial
#' dependence of the measures singles one reading out.
#'
#' @param data A \code{data.frame} (already loaded; for header-aware loading see
#'   [read_data()]).
#' @param name_bonus Optional named list mapping role keys to character vectors
#'   of case-insensitive keyword regex fragments, e.g.
#'   \code{list(group_var = c("treat", "arm"), outcome_binary = c("death"))}.
#'   Recognized keys: \code{group_var}, \code{outcome_continuous},
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
#'       \code{components} score breakdown. Meaning of \code{score}/\code{pct} per role:
#'       \code{group_var}: \code{score = pct = 100 * B} (normalized-entropy ranking;
#'       \code{max_score = 100}; there is no description-length test).
#'       Paired/agreement pairs, ordinal scale items and the outcomes: \code{score}
#'       is in bits, \code{max_score} is \code{NA}, and \code{pct} is the display map
#'       \code{100 * (1 - 2^(-bits/8))}. Survival, wide repeated measures and binary
#'       scale items (ported point scorers): \code{score} and \code{max_score} are
#'       points and \code{pct} is the points ratio.}
#'     \item{value_types}{named character vector of per-column value-type labels.}
#'     \item{potential_pairs}{list of candidate continuous column pairs with
#'       MDL pair bits and shift bits. The paired/agreement split is the sign of the
#'       shift bits, a BIC-type test of a mean difference between the two columns:
#'       \code{agreement_pairs} means "no detectable systematic shift", so a
#'       method-comparison pair with a real calibration bias is labelled
#'       \code{paired_pairs}, and the label depends on n. Exactly one pair is
#'       reported across the two roles.}
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
  if (is.null(names(data)) || anyDuplicated(names(data))) {
    # Robustness (AUREX): duplicate (or absent) column names make `data[[name]]`
    # ambiguous and break value-based addressing; uniquify so every column is
    # addressable. Detection is name-blind, so this never changes role
    # assignment for the normal (unique-name) case (make.unique is a no-op).
    names(data) <- make.unique(if (is.null(names(data))) rep("", ncol(data)) else names(data))
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
